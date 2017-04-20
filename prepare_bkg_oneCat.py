#python g1_exo_doFit_class.py -b -c mu > test.log
#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser
import CMS_lumi, tdrstyle
from array import array


from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooNLLVar, RooAddition, RooProduct, RooConstraintSum, RooCustomizer, RooMinuit, RooArgSet, RooAbsData, RooAbsPdf, RooAbsReal, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite


############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('-b', action='store_true', dest='noX', default=True, help='no X11 windows')
parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--hi', action='store', dest='mlvj_hi', type='float', default=3500, help='dont change atm!')
parser.add_option('--lo', action='store', dest='mlvj_lo', type='float', default=900, help='set lower cut on MWV, mat cause problems')
parser.add_option('-r','--readtrees', action='store_true', dest='read_trees', default=False, help='read data and MC from TTrees, has to be done when range or binning is changed -> takes much longer')
parser.add_option('--noplots', action='store_true', dest='noplots', default=False, help='dont make any plots')


(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/hyperg_2F1_c.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf, RooErfExpDecoPdf



class doFit_wj_and_wlvj:

    def __init__(self, in_channel, in_mj_min=40, in_mj_max=150, in_mlvj_min=900., in_mlvj_max=3500., fit_model="ExpN", fit_model_alter="ExpTail"):

        tdrstyle.setTDRStyle()
        TGaxis.SetMaxDigits(3)
        
        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        ###to save all fitresults
        self.fitresultsmj        = []
        self.fitresultsmlvj        = []
        self.fitresultsfinal        = []
        ### set the channel type --> electron or muon
        self.channel=in_channel;
        self.leg = TLegend(); 
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;
                
        print "########################################################################################"
        print "######## define class: binning, variables, cuts, files and nuissance parameters ########"
        print "########################################################################################"

        ### Set the mj binning for plots
        self.BinWidth_mj=5.;

        ### Set the binning for mlvj plots as a function of the model
        self.BinWidth_mlvj=100.;
            
        ## correct the binning of mj 
        nbins_mj=int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max=in_mj_min+nbins_mj*self.BinWidth_mj;
                   
        ## correct the binning of mlvj 
        nbins_mlvj=int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        in_mlvj_max=in_mlvj_min+nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
        varname = "M_{pruned}"
        rrv_mass_j = RooRealVar("rrv_mass_j",varname,(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV");
        rrv_mass_j.setBins(nbins_mj);

        ## define invariant mass WW variable
        rrv_mass_lvj= RooRealVar("rrv_mass_lvj","M_{WV}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV");
        rrv_mass_lvj.setBins(nbins_mlvj);

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;

        ## create the workspace
        self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_");

        getattr(self.workspace4fit_,"import")(rrv_mass_j);
        getattr(self.workspace4fit_,"import")(rrv_mass_lvj);

        #prepare workspace for unbin-Limit -> just fo the stuff on which running the limit 
        self.workspace4limit_ = RooWorkspace("workspace4limit_","workspace4limit_");

        #define sidebands
        self.mj_sideband_lo_min = int(in_mj_min);
        self.mj_sideband_lo_max = 65
        self.mj_sideband_hi_min = 105
        self.mj_sideband_hi_max = int(in_mj_max);
        self.mj_signal_min = 65
        self.mj_signal_max = 105

        ## zone definition in the jet mass 
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("sig",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

        ## one zone for MWW
        rrv_mass_lvj.setRange("total_region",in_mlvj_min,in_mlvj_max)
        rrv_mass_lvj.setRange('over3500',3500,5000)
        rrv_mass_lvj.setRange('sig',900,3500)

        #prepare the data and mc files --> set the working directory and the files name
        self.file_Directory="InputTrees/"+self.channel+"/";
                 
        #prepare background data and signal samples            
        self.file_data                 = ("tree_data_%s.root"%self.channel);
        self.file_WJets0_mc          = ("tree_WJets_%s.root"%self.channel);
        self.file_WW_mc              = ("tree_WW_%s.root"%self.channel);# WW
        self.file_WZ_mc              = ("tree_WZ_%s.root"%self.channel);# WZ
        self.file_TTbar_mc           = ("tree_TTbar_%s.root"%self.channel);
        self.file_STop_mc            = ("tree_STop_%s.root"%self.channel);

        
        #self.mean_shift = -0.8
        #self.sigma_scale=1.086
        self.mean_shift = -1.294
        self.sigma_scale=0.958

        self.wtagger_label        = 'HPV'
        
        if options.mlvj_hi!=3500:
                raw_input('cuts different from M_WV<3500 not supported atm!')

        self.plotsDir = 'plots_%s_%s_%s_%s' %(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi)
        if not os.path.isdir("cards_%s_%s_%s_%s"%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi)):
                os.system("mkdir cards_%s_%s_%s_%s"%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi));
        self.rlt_DIR="cards_%s_%s_%s_%s/"%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi)

        #result files: The event number, parameters and error write into a txt file. The dataset and pdfs write into a root file


        ## extra text file
        self.file_rlt_txt = self.rlt_DIR+"other_wwlvj_%s_%s.txt"%(self.channel,self.wtagger_label)
        ## workspace for limit
        self.file_rlt_root = self.rlt_DIR+"wwlvj_%s_%s_workspace.root"%(self.channel,self.wtagger_label)
        ## datacard for the ubninned limit
        self.file_datacard_unbin = self.rlt_DIR+"wwlvj_%s_%s_unbin.txt"%(self.channel,self.wtagger_label)
        ## workspace for the binned limit
        self.file_datacard_counting = self.rlt_DIR+"wwlvj_%s_%s_counting.txt"%(self.channel,self.wtagger_label)
        
        self.file_out=open(self.file_rlt_txt,"w");
        self.file_out.write("Welcome:\n");
        self.file_out.close()
        self.file_out=open(self.file_rlt_txt,"a+");

        ## color palette 
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : kGreen+1,
            'WW' : kRed,
            'WZ' : kCyan,
            'STop' : kBlue,
            'TTbar' : kOrange,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        # parameters of data-driven method to get the WJets background event number.
        self.number_WJets_insideband=-1;
        self.datadriven_alpha_WJets_unbin=-1;
        self.datadriven_alpha_WJets_counting=-1;

        #### Set systematic on the Wjets shape   and TTbar due to PS, fitting function etc..
        self.shape_para_error_WJets0 = 1.4;
        self.shape_para_error_alpha  = 1.4;
        self.shape_para_error_TTbar = 2.0;
        self.shape_para_error_VV    = 1.;
        self.shape_para_error_STop  = 1.;
                                                                  
        # shape parameter uncertainty
        self.FloatingParams=RooArgList("floatpara_list");

    #################################################################################################
    #################################################################################################
      
    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1., isalpha=False, ismj=False):
        print "############### draw the legend ########################"
        if left==-1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetLineColor(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.37+x_offset_low, 0.50+y_offset_low, 0.72+x_offset_high, 0.82+y_offset_high, "", "NDC");            
            theLeg.SetName("theLegend");
            if ismj: theLeg = TLegend(0.3715365+x_offset_low,0.505+y_offset_low,0.8526448+x_offset_high,0.845+y_offset_high, "", "NDC"); 
            if TwoCoulum :
                theLeg.SetNColumns(2);
            if isalpha: theLeg = TLegend(0.3944724+x_offset_low,0.4370629+y_offset_low,0.7650754+x_offset_high,0.8374126+y_offset_high, "", "NDC");  
            
        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.05);
        theLeg.SetTextFont(42);

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        if   self.channel == 'el': legHeader="e#nu";
        elif self.channel == 'mu': legHeader="#mu#nu";
        
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          #if objName.find("TPave") != -1: continue
          if objName == "errorband" : objName = "Uncertainty";
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);  objName_before=objName;
            elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "Data, W#rightarrow"+legHeader,"PE");  objName_before=objName;                 
            else: objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        
                   
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="WJets" : objName_before=objName;                 
            else:  objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            if objName.find("TPave") != -1: continue
            if objName == "errorband" : objName = "Uncertainty";
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data()=="STop" : theLeg.AddEntry(theObj, "Single Top","F");
                    #elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ","F");
                    elif TString(objName).Data()=="data" :  objName_before=objName; entryCnt = entryCnt+1; continue ;
                    elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F"); entryCnt = entryCnt+1; continue;
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");

                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
           theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;

    #################################################################################################
    #################################################################################################

    def get_canvas(self,cname,isalpha=False):

       #tdrstyle.setTDRStyle()
       CMS_lumi.lumi_13TeV = "2.3 fb^{-1}"
       CMS_lumi.writeExtraText = 1
       if isalpha:
                       CMS_lumi.extraText = "Simulation\n Preliminary"
       else:
                       CMS_lumi.extraText = "Preliminary"

       iPos = 11
       if( iPos==0 ): CMS_lumi.relPosX = 0.15

       H_ref = 600; 
       W_ref = 800; 
       W = W_ref
       H  = H_ref

       T = 0.08*H_ref
       B = 0.12*H_ref 
       L = 0.12*W_ref
       R = 0.06*W_ref

       canvas = ROOT.TCanvas(cname,cname,W,H)
       canvas.SetFillColor(0)
       canvas.SetBorderMode(0)
       canvas.SetFrameFillStyle(0)
       canvas.SetFrameBorderMode(0)
       canvas.SetLeftMargin( L/W )
       canvas.SetRightMargin( R/W )
       canvas.SetTopMargin( T/H )
       canvas.SetBottomMargin( B/H+0.03 )
       canvas.SetTickx()
       canvas.SetTicky()
       if isalpha:
        canvas.SetTicky(0)
       
       return canvas

    #################################################################################################
    #################################################################################################

    #### just drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0, frompull=0, isalpha=0, fix_axis=0, force_plots=0):
      
        if options.noplots and not force_plots:
          return 0

        print "############### draw the canvas without pull ########################"
        cMassFit = self.get_canvas(in_obj.GetName())#TCanvas("cMassFit","cMassFit", 600,600);

        if fix_axis == 0:
                if frompull and logy :
                    in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/200)
                elif not frompull and logy :
                    in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())

        if is_range:
            h2=TH2D("h2","",100,400,1400,4,0.00001,4);
            h2.Draw();
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045);
        in_obj.GetXaxis().SetTitleOffset(1.15);
        in_obj.GetXaxis().SetLabelSize(0.04);

        in_obj.GetYaxis().SetTitleSize(0.04);
        in_obj.GetYaxis().SetTitleOffset(1.40);
        in_obj.GetYaxis().SetLabelSize(0.04);

        self.leg.SetTextSize(0.031); 
        if isalpha: self.leg.SetTextSize(0.038)

        cMassFit.Update()
        cMassFit.cd()
        if isalpha:
                CMS_lumi.CMS_lumi(cMassFit, 0, 11, 0.075)
        else:
                CMS_lumi.CMS_lumi(cMassFit, 4, 11)
        cMassFit.cd()
        cMassFit.Update()
        cMassFit.RedrawAxis()
        frame = cMassFit.GetFrame()
        frame.Draw()   
        cMassFit.cd()
        cMassFit.Update()
                
        Directory=TString(in_directory); #+"_%02d_%02d/"%(options.cprime,options.BRnew));
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
        else:
            rlt_file.ReplaceAll(".root","");
            rlt_file = rlt_file.Append(".png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".pdf",".root");
        #cMassFit.SaveAs(rlt_file.Data());

        if logy:
            if not isalpha and fix_axis == 0:
                in_obj.GetYaxis().SetRangeUser(1e-3,in_obj.GetMaximum()*200);
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            #cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

    #################################################################################################
    #################################################################################################

    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, rrv_x, datahist, mplot, mplot_pull,ndof,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0,ismj=0, fix_axis=0, force_plots=0):# mplot + pull

        if options.noplots and not force_plots:
          return 0

        print "############### draw the canvas with pull ########################" 
        #@#mplot.GetXaxis().SetTitle(rrv_x.GetTitle() + " (GeV)")
        mplot.GetXaxis().SetTitle("")
        mplot.GetYaxis().SetTitleSize(0.07)
        mplot.GetYaxis().SetTitleOffset(0.9)
        mplot.GetYaxis().SetLabelSize(0.06)
        mplot.GetXaxis().SetLabelSize(0);
        
        cMassFit = self.get_canvas(mplot.GetName())
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
            pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
            pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
            pad3=TPad("pad3","pad3",0.8,0.,1,1);
            pad1.Draw();
            pad2.Draw();
            pad3.Draw();
        else:
            pad1=TPad("pad1","pad1",0.,0. ,1,0.30); #pad1 - pull
            pad2=TPad("pad2","pad2",0.,0.3,1.,1. ); #pad0
            pad2.SetRightMargin(0.1);
            pad2.SetTopMargin(0.1);
            pad2.SetBottomMargin(0.0001);
            pad1.SetRightMargin(0.1)
            pad1.SetTopMargin(0)
            pad1.SetBottomMargin(0.4)   
            pad1.Draw();
            pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
        
        mplot.Draw();

        pad1.cd();
        mplot_pull.Draw("AP");

        medianLine = TLine(mplot.GetXaxis().GetXmin(),0.,mplot.GetXaxis().GetXmax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
        medianLine.Draw()
        mplot_pull.Draw("Psame");
        
        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

        cMassFit.Update()
        pad2.cd()
        CMS_lumi.CMS_lumi(pad2, 4, 11)        
        pad2.cd()
        pad2.Update()
        pad2.RedrawAxis()
        frame = pad2.GetFrame()
        frame.Draw()   
        cMassFit.cd()
        cMassFit.Update()
                        
        ## create the directory where store the plots
        Directory = TString(in_directory);
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
        else:
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","");
            rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());
        
        rlt_file.ReplaceAll(".pdf",".root");
        #cMassFit.SaveAs(rlt_file.Data());

        string_file_name = TString(in_file_name);
        if string_file_name.EndsWith(".root"):
            string_file_name.ReplaceAll(".root","_"+in_model_name);
        else:
            string_file_name.ReplaceAll(".root","");
            string_file_name.Append("_"+in_model_name);

        if logy:
            if fix_axis == 0:
                    mplot.GetYaxis().SetRangeUser(0.002,mplot.GetMaximum()*200);
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            #cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        #self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1,0,fix_axis);

    #################################################################################################
    #################################################################################################

    def calculate_chi2(self,hist,rrv_x,mplot_orig,ndof,ismj):

        pulls = array('d',[])
        print "############### calculate chi2 (new) ########################"
        hpull = mplot_orig.pullHist();
        bins = 0
        bins_ = 0
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
          hist.get(bins_)
          hist.weightError(RooAbsData.SumW2)
          print x,y,bins_,hist.get(bins_).getRealValue(rrv_x.GetName()),hist.weight(),hist.weightError(RooAbsData.SumW2)
          if hist.weight() != 0: pulls.append(y)
          bins_+=1
          
        chi2 = 0
        for p in pulls:
         chi2+=(p*p)
         
        print "Chi2/ndof = %f/%f = %f" %(chi2,ndof,chi2/ndof)
        return chi2,ndof
               
    #################################################################################################
    #################################################################################################

    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist();
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
            hpull.GetPoint(ipoint,x,y);
            #print x,y
            if(y == 0):
                hpull.SetPoint(ipoint,x,10)

        self.gt = ROOT.TH1F(mplot_orig.GetName()+'_pull',mplot_orig.GetName()+'_pull',rrv_x.getBins(),rrv_x.getMin(),rrv_x.getMax());
        self.gt.SetMinimum(-3.999);
        self.gt.SetMaximum(3.999);
        self.gt.SetDirectory(0);
        self.gt.SetStats(0);
        self.gt.SetLineStyle(0);
        self.gt.SetMarkerStyle(20);
        self.gt.GetXaxis().SetTitle(rrv_x.GetTitle() + " (GeV)");
        self.gt.GetXaxis().SetLabelFont(42);
        self.gt.GetXaxis().SetLabelOffset(0.02);
        self.gt.GetXaxis().SetLabelSize(0.15);
        self.gt.GetXaxis().SetTitleSize(0.15);
        self.gt.GetXaxis().SetTitleOffset(1.2);
        self.gt.GetXaxis().SetTitleFont(42);
        self.gt.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
        self.gt.GetYaxis().CenterTitle(True);
        self.gt.GetYaxis().SetNdivisions(205);
        self.gt.GetYaxis().SetLabelFont(42);
        self.gt.GetYaxis().SetLabelOffset(0.007);
        self.gt.GetYaxis().SetLabelSize(0.15);
        self.gt.GetYaxis().SetTitleSize(0.15);
        self.gt.GetYaxis().SetTitleOffset(0.35);
        self.gt.GetYaxis().SetTitleFont(42);
        hpull.SetHistogram(self.gt)

        return hpull

    #################################################################################################
    #################################################################################################

    #set uncertainties for data as recommended by the statistics commitee
    def getData_PoissonInterval(self,data_obs,mplot):
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        datahist   = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone");
        data_histo = datahist.createHistogram("histo_data",rrv_x) ;
        data_histo.SetName("data");
        data_plot  = RooHist(data_histo);
        data_plot.SetMarkerStyle(20);
        data_plot.SetMarkerSize(1);
        
        alpha = 1 - 0.6827;
        for iPoint  in range(data_plot.GetN()):
                  N = data_plot.GetY()[iPoint];
                  if N==0 : 
                        L = 0;
                  else : 
                        L = (ROOT.Math.gamma_quantile(alpha/2,N,1.));
                  U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1);
                  data_plot.SetPointEYlow(iPoint, N-L);
                  data_plot.SetPointEYhigh(iPoint,U-N);
                  data_plot.SetPointEXlow(iPoint,0);        
                  data_plot.SetPointEXhigh(iPoint,0);        
        
        mplot.addPlotable(data_plot,"PE");

    #################################################################################################
    #################################################################################################

    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self,label, mlvj_region):
        return self.workspace4fit_.pdf("model"+label+mlvj_region+"_"+self.channel+"_mlvj");
    
    #################################################################################################
    #################################################################################################        
        
    ### get an mj model from the workspace given the label
    def get_mj_Model(self,label):
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")

    #################################################################################################
    #################################################################################################

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region="_sig"):
        print "########### Fixing a general mlvj model  ############"
        rdataset_General_mlvj = self.workspace4fit_.data("rdataset%s%s_%s_mlvj"%(label, mlvj_region,self.channel))
        model_General = self.get_mlvj_Model(label,mlvj_region);
        rdataset_General_mlvj.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_mlvj);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return model_General
        


    #################################################################################################
    #################################################################################################

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_mj_Model(self, label ):
        print "########### Fixing a general mj model  ############"
        rdataset_General_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_General = self.get_mj_Model(label);
        rdataset_General_mj.Print();
        model_General.Print();
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_mj);
        par=parameters_General.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName())
            if ( paraName.Contains("rrv_c_ErfExp_WJets") or paraName.Contains("rrv_p0_User1_WJets")):
                     param.setConstant(kFALSE);
                     param.Print();
            else:
                    param.setConstant(kTRUE);
            param=par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    #################################################################################################
    #################################################################################################

    ### fix only the WJets model --> for extended pdf (just fix shape parameters of width, offset of ErfExp and p1 of User1 function
    def get_WJets_mj_Model(self,label):
        print "########### Fixing only the WJets mj Shape --> just the printed parameters  ############"
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_WJets = self.get_mj_Model(label);
        rdataset_WJets_mj.Print();
        model_WJets.Print();
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mj);
        par=parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            if ( paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or paraName.Contains("rrv_p1_User1_WJets")):
             param.setConstant(kTRUE);
             param.Print();
            else:
             param.setConstant(0);
            param=par.Next()
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    #################################################################################################
    #################################################################################################

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="total_region",mass_spectrum="_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj  ############"        
        rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
        model = self.get_mlvj_Model(label,mlvj_region);
        model.Print();
        rdataset.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param=par.Next()

    #################################################################################################
    #################################################################################################

    ### fix a pdf in a different way --> for RooAbsPdf 
    def fix_Pdf(self,model_pdf,argset_notparameter):
        print "########### Fixing a RooAbsPdf for mlvj or mj  ############"        
        parameters = model_pdf.getParameters(argset_notparameter);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()


    #################################################################################################
    #################################################################################################

    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters          
    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):
        if TString(mass_spectrum).Contains("_mj"): rrv_x = self.workspace4fit_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4fit_.var("rrv_mass_lvj");

        ## ExpN pdf for W+jets bkg fit
        if in_model_name == "ExpN":
            
            print "########### ExpN funtion for W+jets mlvj ############"
            rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-2e-3,-1e-1,-1e-5);
            rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 4e3, -1e4, 1e4);
            if rrv_x.getMin() == 700:
               rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 0, -10000, 10000);               

            model_pdf = ROOT.RooExpNPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ExpN, rrv_n_ExpN);
                                                                                                             
        ## levelled exp for W+jets bkg fit
        if in_model_name == "ExpTail":
            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if self.wtagger_label.find("LP") != -1:
             rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
             rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1e-1,-1.e-2,1e6);
            else:
                if self.channel == "el" :
                 if ismc == 1 and label_tstring.Contains("sb"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 139,0.,355);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2e-2,-1.e-2,5.5e-2);                     
                 elif ismc == 1 and label_tstring.Contains("sig"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 162,18,395);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1.6e-2,-1.e-2,5.5e-2);
                 elif ismc == 0 :  
                     rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,70,240);
                     rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);
                           
                if self.channel == "mu" or self.channel == "em":
                 if ismc == 1 and label_tstring.Contains("sb"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-1e-2,7.5e-2);                                      
                 elif ismc == 1 and label_tstring.Contains("sig"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 110,20,500);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1,7.5e-2);
                 elif ismc == 0 :  
                     rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,40,280);
                     rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);    
      

            rrv_s_ExpTail.Print()     
            rrv_a_ExpTail.Print()     
            model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

        ## sum of two exponential 
        if in_model_name == "Exp" or in_model_name == "Exp_sr":
            print "########### Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,-0.05,-0.1,0.);
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp);

        ## Erf times Exp for mj spectrum
        if in_model_name == "ErfExp" :
            print "########### Erf*Exp for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,15.,-30,30);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,60.,30.,150);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,50.,30, 100.);
            #model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);
            model_pdf         = ROOT.RooErfExpDecoPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## User1 function 
        if in_model_name == "User1":
            print "########### User 1 Pdf  for mlvj fit ############"
            rrv_p0     = RooRealVar("rrv_p0_User1"+label+"_"+self.channel,"rrv_p0_User1"+label+"_"+self.channel, 12, 10, 30);
            if self.wtagger_label=="HP": #change this!
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel,"rrv_p1_User1"+label+"_"+self.channel, -4, -9, -2);
            else:
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel,"rrv_p1_User1"+label+"_"+self.channel, -2.5, -4, 0.);
            model_pdf=RooUser1Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_p0,rrv_p1);

        ## Exp+Gaus or mj spectrum
        if in_model_name == "ExpGaus":
            print "########### Exp + Gaus for mj  fit  ############"
            rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,0.05,-0.2,0.2);
            exp             = ROOT.RooExponential("exp"+label+"_"+self.channel,"exp"+label+"_"+self.channel,rrv_x,rrv_c_Exp);

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,78,90);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,15);
            rrv_high        = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.56,0.,1.);
            gaus            = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            model_pdf       = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high),1)


        ## Erf*Exp + 2Gaus for mj spectrum
        ###TTBARMJFIT
        if in_model_name == "2Gaus_ErfExp":

            print "########### 2Gaus + Erf*Exp for mj fit  ############"
            mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            #frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;
            frac_tmp       = 0; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp, deltamean_tmp-deltamean_tmp_err*10, deltamean_tmp+deltamean_tmp_err*10);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*10, scalesigma_tmp+scalesigma_tmp_err*10);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel,"rrv_frac_2gaus"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            c0_tmp     = -2.9893e-02 ; c0_tmp_err     = 6.83e-03;
            offset_tmp = 7.9350e+01  ; offset_tmp_err = 9.35e+00;
            width_tmp  = 3.3083e+01  ; width_tmp_err  = 2.97e+00;

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-10, width_tmp+10);
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum,"erfexp"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel, 0.5,0.,1.);
            #use only one gauss
            #model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1),RooArgList(rrv_frac),1)

        #2Gaus for WW
        if in_model_name == "2GausWW":
            
            print "########### 2Gaus for mj fit  ############"
            mean1_tmp      = 8.385e+01; mean1_tmp_err      = 1.63e-01;
            #mean1_tmp      = 7.3141e+01; mean1_tmp_err      = 1.63e-01; @@@ JEN            
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-6, mean1_tmp+6);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-6,sigma1_tmp+6 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,0.,-8,15);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*6, scalesigma_tmp+scalesigma_tmp_err*6);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac1         = RooRealVar("rrv_frac1"+label+"_"+self.channel,"rrv_frac1"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*10, frac_tmp+frac_tmp_err*10);
            gausguas_1         = RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum,"gausguas_1"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)
            
            model_pdf         = gausguas_1.clone("model_pdf"+label+"_"+self.channel+mass_spectrum)
        #2Gaus for WZ
        if in_model_name == "2GausWZ":
            
            print "########### 2Gaus for mj fit  ############"
            mean1_tmp      = 9.11876e+01; mean1_tmp_err      = 1.63e-01;
            #mean1_tmp      = 7.3141e+01; mean1_tmp_err      = 1.63e-01; @@@ JEN            
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-6, mean1_tmp+6);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-6,sigma1_tmp+6 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,0.,-8,10);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*6, scalesigma_tmp+scalesigma_tmp_err*6);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac1         = RooRealVar("rrv_frac1"+label+"_"+self.channel,"rrv_frac1"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*10, frac_tmp+frac_tmp_err*10);
            gausguas_1         = RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum,"gausguas_1"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)
            
            model_pdf         = gausguas_1.clone("model_pdf"+label+"_"+self.channel+mass_spectrum)
            
            
        ## return the pdf
        getattr(self.workspace4fit_,"import")(model_pdf)
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)

    #################################################################################################
    #################################################################################################

    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500):

      ##### define an extended pdf from a standard Roofit One
      print " "
      print "###############################################"
      print "## Make model : ",label," ",in_model_name,"##";
      print "###############################################"
      print " "

      ## call the make RooAbsPdf method
      rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum,"rrv_number"+label+"_"+self.channel+mass_spectrum,area_init_value,0.,1e7);
      model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList,ismc_wjet)
      print "######## Model Pdf ########"        
      model_pdf.Print();
      rrv_number.Print();
      
      ## create the extended pdf
      model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
      print "######## Model Extended Pdf ########"        

      #### put all the parameters ant the shape in the workspace
      getattr(self.workspace4fit_,"import")(rrv_number)
      getattr(self.workspace4fit_,"import")(model) 
      self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();
      ## return the total extended pdf
      return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum);

    #################################################################################################
    #################################################################################################

    ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    ###only W+jets is fitted, other bkgs are taken from fits to MC
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0):

        print "############### Fit mlvj in mj sideband: ",label," ",mlvj_region,"  ",mlvj_model," ##################"
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset_data_mlvj = self.workspace4fit_.data("rdataset_data%s_%s_mlvj"%(mlvj_region,self.channel))

        ## get the minor component shapes in the sb low
        model_WW_backgrounds    = self.get_General_mlvj_Model("_WW","_sb").clone("WW")
        number_WW_sb_mlvj            = self.workspace4fit_.var("rrv_number_WW_sb_%s_mlvj"%(self.channel))
        model_WZ_backgrounds    = self.get_General_mlvj_Model("_WZ","_sb").clone("WZ")
        number_WZ_sb_mlvj            = self.workspace4fit_.var("rrv_number_WZ_sb_%s_mlvj"%(self.channel))
        model_TTbar_backgrounds = self.get_General_mlvj_Model("_TTbar","_sb").clone("TTbar")
        number_TTbar_sb_mlvj         = self.workspace4fit_.var("rrv_number_TTbar_sb_%s_mlvj"%(self.channel))
        model_STop_backgrounds  = self.get_General_mlvj_Model("_STop","_sb").clone("STop")
        number_STop_sb_mlvj          = self.workspace4fit_.var("rrv_number_STop_sb_%s_mlvj"%(self.channel))

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_sb_from_fitting"%(label), mlvj_model,"_mlvj");
        model_pdf_WJets.Print();
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_sb = self.workspace4fit_.var("rrv_number%s_sb_%s_mlvj"%(label,self.channel)).clone("rrv_number%s_sb_from_fitting_%s_mlvj"%(label,self.channel));
        model_WJets =RooExtendPdf("model%s_sb_from_fitting_%s_mlvj"%(label,self.channel),"model%s_sb_from_fitting_%s_mlvj"%(label,self.channel),model_pdf_WJets,number_WJets_sb);
        model_pdf_WJets.Print();
        number_WJets_sb.Print()

        ## Add the other bkg component fixed to the total model
        model_data = RooAddPdf("model_data%s%s_mlvj"%(label,mlvj_region),"model_data%s%s_mlvj"%(label,mlvj_region),RooArgList(model_WJets,model_WW_backgrounds,model_WZ_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds));
        
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        self.fitresultsmlvj.append(rfresult)
        getattr(self.workspace4fit_,"import")(model_WJets)

        model_WJets.Print();
        model_WJets.getParameters(rdataset_data_mlvj).Print("v");
        self.workspace4fit_.pdf("model_pdf%s_sb_%s_mlvj"%(label,self.channel)).getParameters(rdataset_data_mlvj).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_sb_mlvj = RooRealVar("rrv_number_data_sb_%s_mlvj"%(self.channel),"rrv_number_data_sb_%s_mlvj"%(self.channel),
                                                 self.workspace4fit_.var("rrv_number_TTbar_sb_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number_STop_sb_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number_WW_sb_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number_WZ_sb_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number%s_sb_from_fitting_%s_mlvj"%(label,self.channel)).getVal() );
        ### anti-correlation between different components not taken into account, but still only hardly visible in plot
        rrv_number_data_sb_mlvj.setError( TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_sb_from_fitting_%s_mlvj"%(label,self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number%s_sb_from_fitting_%s_mlvj"%(label,self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_TTbar_sb_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_TTbar_sb_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_STop_sb_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_STop_sb_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_WW_sb_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_WW_sb_%s_mlvj"%(self.channel)).getError()+                                                        self.workspace4fit_.var("rrv_number_WZ_sb_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_WZ_sb_%s_mlvj"%(self.channel)).getError()));

        getattr(self.workspace4fit_,"import")(rrv_number_data_sb_mlvj)

        ### plot for WJets default + default shape
        if TString(label).Contains("_WJets0"):

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(rrv_mass_lvj.getBins()));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0) );

            model_data.plotOn(mplot, RooFit.Components("model%s_sb_from_fitting_%s_mlvj,TTbar,STop,WW,WZ"%(label,self.channel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(),RooAbsReal.NumEvent)) ;
            model_data.plotOn(mplot, RooFit.Components("TTbar,STop,WW,WZ"),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(),RooAbsReal.NumEvent)) ;
            model_data.plotOn(mplot, RooFit.Components("WW,WZ,STop"), RooFit.Name("WW"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(),RooAbsReal.NumEvent));
            model_data.plotOn(mplot, RooFit.Components("WZ,STop"), RooFit.Name("WZ"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WZ"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(),RooAbsReal.NumEvent)); 
            model_data.plotOn(mplot, RooFit.Components("STop"), RooFit.Name("STop"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(),RooAbsReal.NumEvent));

            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_sb_from_fitting_%s_mlvj,TTbar,STop,WW,WZ"%(label,self.channel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(), RooAbsReal.NumEvent)) ;
            model_data.plotOn(mplot, RooFit.Components("TTbar,STop,WW,WZ"),RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(), RooAbsReal.NumEvent)) ;
            model_data.plotOn(mplot, RooFit.Components("WW,WZ,STop"), RooFit.Name("WW_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(), RooAbsReal.NumEvent));
            model_data.plotOn(mplot, RooFit.Components("WZ,STop"), RooFit.Name("WZ_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(), RooAbsReal.NumEvent));
            model_data.plotOn(mplot, RooFit.Components("STop"), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_mlvj.getVal(), RooAbsReal.NumEvent));
 
            ### draw the error band 
            draw_error_band(rdataset_data_mlvj, model_data,self.workspace4fit_.var("rrv_number_data_sb_%s_mlvj"%(self.channel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible());
            self.getData_PoissonInterval(rdataset_data_mlvj,mplot);

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);  

            ### Add the legend to the plot 
            #self.leg=self.legend4Plot(mplot,0,1,0., 0.16, 0.16, 0.1);
            self.leg=self.legend4Plot(mplot,0,1,0.25,0.,0.1,0.,0);
            mplot.addObject(self.leg)

            ### calculate the chi2
            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            #print mplot.chiSquare();
            #print "#################### nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            ### write the result in the output
            #self.file_out.write("\n fit_mlvj_in_Mj_sideband: nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof ) );

            ### get the pull plot and store the canvas
            
            mplot_pull = self.get_pull(rrv_mass_lvj,mplot);
            parameters_list = model_data.getParameters(rdataset_data_mlvj);
            
            mplot.GetXaxis().SetTitle(rrv_mass_lvj.GetTitle() + " (GeV)")
            mplot.GetYaxis().SetRangeUser(2e-3,2e5)
            #@#self.draw_canvas_with_pull( rrv_mass_lvj,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), "m_lvj_sb%s"%(label),"",1,1)
            self.draw_canvas_with_pull( rrv_mass_lvj,rdataset_data_mlvj,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), "m_lvj_sb%s"%(label),"",1,1,0,1)
            

        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sb_from_fitting_mlvj"%(label));
        purity = self.wtagger_label[0]+self.wtagger_label[1]
        Deco      = PdfDiagonalizer("Deco%s_sb_from_fitting_%s_%s_mlvj_13TeV"%(label,self.channel,purity),wsfit_tmp,rfresult);
        print"#################### diagonalize data sideband fit "
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets);
        print "##################### workspace for decorrelation ";
        wsfit_tmp.Print("v");
        print "##################### original  parameters ";
        model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
        print "##################### original  decorrelated parameters ";
        model_pdf_WJets_deco.getParameters(rdataset_data_mlvj).Print("v");
        print "##################### original  pdf ";
        model_pdf_WJets.Print();
        print "##################### decorrelated pdf ";
        model_pdf_WJets_deco.Print();

        getattr(self.workspace4fit_,"import")(model_pdf_WJets_deco);

        #### Call the alpha evaluation in automatic
        self.get_WJets_mlvj_correction_sb_to_sig(label,mlvj_model);

        ### Fix the pdf of signal, TTbar, STop WW and WZ in the signal region 
        self.fix_Model("_TTbar","_sig","_mlvj")
        self.fix_Model("_STop","_sig","_mlvj")
        self.fix_Model("_WW","_sig","_mlvj")
        self.fix_Model("_WZ","_sig","_mlvj")

    #################################################################################################
    #################################################################################################

    ##### Function that calculates the normalization inside the mlvj signal region
    def get_mlvj_normalization_insignalregion(self, label, model_name=""):
        
        print "############### get mlvj normalization inside SR ",label," ",model_name," ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_sig"+"_"+self.channel+"_mlvj");
        else:
            model = self.workspace4fit_.pdf(model_name);

        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj");
        fullInt   = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj) );
        signalInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("sig"));
        highMassInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("high_mass"));
        over3500Int = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("over3500"))
        
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        highMassInt_val = highMassInt.getVal()/fullInt_val 
        over3500Int_val = over3500Int.getVal()/fullInt_val

        ## integral in the signal region
        print "######### integral in SR: ",label+"signalInt=%s"%(signalInt_val)

        print "####### Events Number in MC Dataset:"
        self.workspace4fit_.var("rrv_number_dataset_sig"+label+"_"+self.channel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").Print();
        if 'WJets1' not in label:
                print "########## Events Number get from fit:"
                rrv_tmp                = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_postfit_sig");
                rrv_tmp_pre        = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_prefit_sig");

                print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)
                rrv_tmp.setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_postfit_sig").getVal()*(signalInt_val-over3500Int_val))
                rrv_tmp_pre.setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_prefit_sig").getVal()*(signalInt_val-over3500Int_val))
        else:
                rrv_tmp                = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj")

        #### store the info in the output file
        self.file_out.write( "\n%s++++++++++++++++++get_mlvj_normalization_insignalregion++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in All Region from dataset : %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset: %s"%(self.workspace4fit_.var("rrv_number_dataset_sig"+label+"_"+self.channel+"_mlvj").getVal()) )
        self.file_out.write( "\nRatio sig/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_sig"+label+"_"+self.channel+"_mlvj").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal() ) )
        self.file_out.write( "\nEvents Number in All Region from fitting : %s\n"%(rrv_tmp.getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting: %s\n"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nEvents Number in High Mass Region from fitting: %s\n"%(rrv_tmp.getVal()*highMassInt_val) )
        self.file_out.write( "\nRatio sig/all_range from fitting :%s"%(signalInt_val ) )
        if not self.workspace4fit_.var("rrv_number_fitting_sig"+label+"_"+self.channel+"_mlvj"):
            rrv_number_fitting_sig_mlvj = RooRealVar("rrv_number_fitting_sig"+label+"_"+self.channel+"_mlvj","rrv_number_fitting_sig"+label+"_"+self.channel+"_mlvj", rrv_tmp.getVal()*signalInt_val );
            getattr(self.workspace4fit_,"import")(rrv_number_fitting_sig_mlvj);
        else :
            self.workspace4fit_.var("rrv_number_fitting_sig"+label+"_"+self.channel+"_mlvj").setVal(rrv_tmp.getVal()*signalInt_val);

        self.workspace4fit_.var("rrv_number_fitting_sig"+label+"_"+self.channel+"_mlvj").Print();

    #################################################################################################
    #################################################################################################

    ### method to get the alpha function to extrapolate the wjets in the signal region
    def get_WJets_mlvj_correction_sb_to_sig(self,label, mlvj_model):

        print" ############# get the extrapolation function alpha from MC : ",label,"   ",mlvj_model," ###############";          

        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here 
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        rdataset_WJets_sb_mlvj          = self.workspace4fit_.data("rdataset4fit%s_sb_%s_mlvj"%(label,self.channel))
        rdataset_WJets_sig_mlvj         = self.workspace4fit_.data("rdataset4fit%s_sig_%s_mlvj"%(label,self.channel))

        ### create a frame for the next plots 
        mplot = rrv_x.frame(RooFit.Title("correlation_pdf"), RooFit.Bins(rrv_x.getBins())) ;
        mplot.GetYaxis().SetTitle("F_{W+jets}^{SR,MC},F_{W+jets}^{SB,MC} (Arbitrary units)");
                
        ### define alpha function depending on used signal model
        if mlvj_model=="ExpN":
            rrv_c_sb  = self.workspace4fit_.var("rrv_c_ExpN%s_sb_%s"%(label,self.channel));
            rrv_n_sb  = self.workspace4fit_.var("rrv_n_ExpN%s_sb_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_ExpN%s_%s"%(label,self.channel),"rrv_delta_c_ExpN%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_sig_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal(),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_sig_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_sig_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )
            rrv_delta_n = RooRealVar("rrv_delta_n_ExpN%s_%s"%(label,self.channel),"rrv_delta_n_ExpN%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_sig_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal(),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_sig_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal()-4*rrv_n_sb.getError(),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_sig_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal()+4*rrv_n_sb.getError() )
            correct_factor_pdf = RooExpNPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c, rrv_delta_n);
 
        if mlvj_model=="ExpTail":
            rrv_s_sb =self.workspace4fit_.var("rrv_s_ExpTail%s_sb_%s"%(label,self.channel));
            rrv_a_sb =self.workspace4fit_.var("rrv_a_ExpTail%s_sb_%s"%(label,self.channel));

            rrv_delta_s = RooRealVar("rrv_delta_s_ExpTail%s_%s"%(label,self.channel),"rrv_delta_s_ExpTail%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_sig_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal(),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_sig_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal()-4*rrv_s_sb.getError(),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_sig_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal()+4*rrv_s_sb.getError() )
            rrv_delta_a = RooRealVar("rrv_delta_a_ExpTail%s_%s"%(label,self.channel),"rrv_delta_a_ExpTail%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_sig_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal(),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_sig_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal()-4*rrv_a_sb.getError(),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_sig_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal()+4*rrv_a_sb.getError() )
                     
            rrv_a_sr = RooFormulaVar("rrv_a_ExpTail_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_a_sb, rrv_delta_a ) );
            rrv_s_sr = RooFormulaVar("rrv_s_ExpTail_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_s_sb, rrv_delta_s ) );
            correct_factor_pdf = RooAlpha4ExpTailPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_s_sr, rrv_a_sr, rrv_s_sb, rrv_a_sb);

        ### define the category and do the simultaneous fit taking the combined dataset of events in mlvj sideband (data) and signal region (MC)
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("sig");
        combData4fit = self.workspace4fit_.data("combData4fit%s_%s"%(label,self.channel));

        model_pdf_sb_WJets      = self.workspace4fit_.pdf("model_pdf%s_sb_%s_mlvj"%(label,self.channel));
        model_pdf_sig_WJets         = RooProdPdf("model_pdf%s_sig_%s_mlvj"%(label,self.channel),"model_pdf%s_sig_%s_mlvj"%(label,self.channel) ,model_pdf_sb_WJets,correct_factor_pdf);

        simPdf = RooSimultaneous("simPdf","simPdf",data_category);
        simPdf.addPdf(model_pdf_sb_WJets,"sideband");
        simPdf.addPdf(model_pdf_sig_WJets,"sig");
        simPdf.fitTo(combData4fit, RooFit.SumW2Error(kTRUE));
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
        #self.fitresultsfinal.append(rfresult)

        ### Decorrelate the parameters in the alpha shape
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sim_mlvj"%(label));
        print "############### diagonalizer alpha ";
        Deco      = PdfDiagonalizer("Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label),wsfit_tmp,rfresult);
        correct_factor_pdf_deco = Deco.diagonalize(correct_factor_pdf);
        print "##################### workspace for decorrelation ";
        wsfit_tmp.Print("v");
        print "##################### original  parameters ";
        correct_factor_pdf.getParameters(rdataset_WJets_sig_mlvj).Print("v");
        print "##################### original  decorrelated parameters ";
        correct_factor_pdf_deco.getParameters(rdataset_WJets_sig_mlvj).Print("v");
        print "##################### original  pdf ";
        correct_factor_pdf.Print();
        print "##################### decorrelated pdf ";
        correct_factor_pdf_deco.Print();

        getattr(self.workspace4fit_,"import")(correct_factor_pdf_deco);
                     

        ### Total plot shape in sb, sr and alpha
        model_pdf_sb_WJets.plotOn(mplot,RooFit.Name("Sideband"),RooFit.LineStyle(kDashed));
        model_pdf_sig_WJets.plotOn(mplot, RooFit.LineColor(kRed), RooFit.Name('Signal Region'));
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("Transfer function #alpha"));

        ### plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
        if TString(label).Contains("_WJets0"):
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha: Alternate PS") );

            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange+7), RooFit.LineStyle(7),RooFit.Name("#alpha: Alternate Function") );

        paras=RooArgList();

        if mlvj_model=="ExpN" or mlvj_model=="ExpTail":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
        
        if TString(label).Contains("_WJets0") or TString(label).Contains("_WJets1"): ### draw error band at 1 and 2 sigma using the decorrelated shape
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGreen,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,2 ,mplot,kYellow,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGreen,"F",3001,"#alpha_invisible #pm",20,400);
            
        ### plot on the same canvas
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha_invisible"))

        if TString(label).Contains("_WJets0") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange+7), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function"));

        elif TString(label).Contains("_WJets1") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange-3), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        ### Add the legend
        self.leg= self.legend4Plot(mplot,1,0, 0, 0., 0., -0.1, 0., True);
        mplot.addObject(self.leg);
        
        ## set the Y axis in arbitrary unit 
        tmp_y_max=0.125
        tmp_y_min=1e-5
        mplot.GetYaxis().SetRangeUser(tmp_y_min,tmp_y_max);

        #### Draw another axis with the real value of alpha
        model_pdf_sb_WJets.getVal(RooArgSet(rrv_x))
        model_pdf_sig_WJets.getVal(RooArgSet(rrv_x))
        correct_factor_pdf_deco.getVal(RooArgSet(rrv_x))
        tmp_alpha_ratio         = ( model_pdf_sig_WJets.getVal(RooArgSet(rrv_x))/model_pdf_sb_WJets.getVal(RooArgSet(rrv_x)) );
        tmp_alpha_pdf           = correct_factor_pdf_deco.getVal(RooArgSet(rrv_x)) * mplot.getFitRangeBinW(); ## value of the pdf in each point
        tmp_alpha_scale         = tmp_alpha_ratio/tmp_alpha_pdf;

        #add alpha scale axis
        axis_alpha=TGaxis( rrv_x.getMax(), 0, rrv_x.getMax(), tmp_y_max, tmp_y_min * tmp_alpha_scale, tmp_y_max*tmp_alpha_scale, 510, "+L" ); #-,-+,+,L
        axis_alpha.SetTitle("Transfer function #alpha");
        axis_alpha.SetTitleOffset(0.75);
        axis_alpha.SetTitleSize(0.05);
        axis_alpha.SetLabelSize(0.045);
        axis_alpha.SetTitleFont(42);
        axis_alpha.SetLabelFont(42);
        #axis_alpha.RotateTitle(1);
        mplot.addObject(axis_alpha);

        self.draw_canvas(mplot,"%s/other/"%(self.plotsDir),"correction_pdf%s_%s_M_lvj_sig_to_sideband"%(label,mlvj_model),0,1,0,1);

        #@#make the same plot with log-scale
        if 'WJets0' in label:
            mplot2 = rrv_x.frame(RooFit.Title("correlation_pdf_log"), RooFit.Bins(rrv_x.getBins())) ;
            mplot3 = rrv_x.frame(RooFit.Title("correlation_pdf_alpha"), RooFit.Bins(rrv_x.getBins()))
            mplot2.GetYaxis().SetTitle("F_{W+jets}^{SR,MC},F_{W+jets}^{SB,MC} (Arbitrary units)");
            model_pdf_sig_WJets.plotOn(mplot2, RooFit.LineColor(kRed), RooFit.Name("Signal Region"));
            model_pdf_sb_WJets.plotOn(mplot2,RooFit.Name("Sideband"),RooFit.LineStyle(kDashed));
           
            correct_factor_pdf_deco.plotOn(mplot3, RooFit.LineColor(kBlack),RooFit.Name("Transfer function #alpha"));
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot3,kGreen,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,2 ,mplot3,kYellow,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot3,kGreen,"F",3001,"#alpha_invisible #pm",20,400);
            correct_factor_pdf_deco.plotOn(mplot3, RooFit.LineColor(kBlack),RooFit.Name("Transfer function #alpha"));
             
            c1        = TCanvas('alpha_log_plot','aplha_log_plot',800,600)
            c1.SetRightMargin(0.1)
            leg        = self.legend4Plot(mplot,1,0, 0, 0.1, 0., 0.075, 0., True)
            pad         = TPad('pad','pad',0,0,1,1)
            pad.SetRightMargin(0.1)
            pad.SetFillColor(0)
            pad.SetLogy()
            pad.SetTicky(0)
            mplot2.GetYaxis().SetRangeUser(1e-5,10)
            pad.Draw()
            pad.cd()
            mplot2.addObject(leg)
            mplot2.Draw()
            
            c1.cd()
            tmp=''
            pad_log = TPad('log','log',0,0,1,1)
            pad_log.SetFillStyle(4000)
            pad_log.SetFillColor(0)
            pad_log.SetFrameFillStyle(4000)
            pad_log.SetRightMargin(0.1)
            pad_log.Draw()
            pad_log.cd()
            mplot3.addObject(axis_alpha)
            tmp=''
            mplot3.GetYaxis().SetLabelOffset(99)
            mplot3.GetYaxis().SetTitleOffset(99)
            mplot3.GetYaxis().SetRangeUser(tmp_y_min,tmp_y_max)
            mplot3.GetYaxis().SetNdivisions(0)
            mplot3.Draw()
            c1.SaveAs(self.plotsDir+'/other/alpha_%s_%s.png'%(self.channel,self.wtagger_label))
            c1.SaveAs(self.plotsDir+'/other/alpha_%s_%s.pdf'%(self.channel,self.wtagger_label))
            
            pad.Delete()
            pad_log.Delete()
            #@#

        correct_factor_pdf_deco.getParameters(rdataset_WJets_sb_mlvj).Print("v");
        model_pdf_WJets_sb_from_fitting_mlvj_Deco = self.workspace4fit_.pdf("model_pdf%s_sb_from_fitting_%s_mlvj_Deco%s_sb_from_fitting_%s_%s_mlvj_13TeV"%(label,self.channel,label, self.channel,self.wtagger_label[0]+self.wtagger_label[1]));
        model_pdf_WJets_sb_from_fitting_mlvj_Deco.Print("v");
        model_pdf_WJets_sb_from_fitting_mlvj = self.workspace4fit_.pdf("model_pdf%s_sb_from_fitting_%s_mlvj"%(label,self.channel))
        model_pdf_WJets_sb_from_fitting_mlvj.Print("v");

        ### Wjets shape in the SR correctedfunction * sb 
        model_pdf_WJets_sig_after_correct_mlvj = RooProdPdf("model_pdf%s_sig_%s_after_correct_mlvj"%(label,self.channel),"model_pdf%s_sig_%s_after_correct_mlvj"%(label,self.channel),model_pdf_WJets_sb_from_fitting_mlvj_Deco,self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label)) );
        model_pdf_WJets_sig_after_correct_mlvj.Print("v")
        #@#also with undeccorrelated pdf
        model_pdf_WJets_sig_undeco_mlvj = RooProdPdf("model_pdf%s_sig_%s_undeco_mlvj"%(label,self.channel),"model_pdf%s_sig_%s_undeco_mlvj"%(label,self.channel),model_pdf_sb_WJets,self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label)) );
        model_pdf_WJets_sig_undeco_mlvj.Print("V")
        ### fix the parmaters and import in the workspace
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_sig_after_correct_mlvj)
        getattr(self.workspace4fit_,'import')(model_pdf_WJets_sig_undeco_mlvj,RooFit.RecycleConflictNodes())

        ##### calculate the normalization and alpha for limit datacard
        self.workspace4fit_.var("rrv_number%s_sig_%s_mlvj"%(label,self.channel)).Print();
        self.workspace4fit_.var("rrv_number%s_sig_%s_mlvj"%(label,self.channel)).setConstant(kTRUE);

    #################################################################################################
    #################################################################################################

    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0): # to get the normalization of WJets in sig

        self.mj_prefit_plot()

    #################################################################################################
    #################################################################################################

    def mj_prefit_plot(self): 

        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rdataset_data_mj=self.workspace4fit_.data("rdataset_data_%s_mj"%(self.channel))

        ### Fix TTbar, VV and STop
        model_TTbar = self.get_General_mj_Model("_TTbar");
        model_STop  = self.get_General_mj_Model("_STop");
        model_WW    = self.get_General_mj_Model("_WW");
        model_WZ    = self.get_General_mj_Model("_WZ");
        model_WJets = self.get_General_mj_Model('_WJets0');

        ## Total Pdf and fit only in sideband 
        model_data = RooAddPdf("model_data_%s_mj"%(self.channel),"model_data_%s_mj"%(self.channel),RooArgList(model_WJets,model_WW,model_WZ,model_TTbar,model_STop))
        getattr(self.workspace4fit_,"import")(model_data)

        ## Total numver of event 
        rrv_number_data_mj = RooRealVar("rrv_number_data_%s_mj"%self.channel,"rrv_number_data_%s_mj"%self.channel,
                                         self.workspace4fit_.var("rrv_number_TTbar_%s_mj"%self.channel).getVal()+
                                         self.workspace4fit_.var("rrv_number_STop_%s_mj"%self.channel).getVal()+
                                         self.workspace4fit_.var("rrv_number_WW_%s_mj"%self.channel).getVal()+
                                         self.workspace4fit_.var("rrv_number_WZ_%s_mj"%self.channel).getVal()+
                                         self.workspace4fit_.var("rrv_number_WJets0_%s_mj"%self.channel).getVal());

        rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number_TTbar_%s_mj"%self.channel).getError()**2+
                                               self.workspace4fit_.var("rrv_number_STop_%s_mj"%self.channel).getError()**2+
                                               self.workspace4fit_.var("rrv_number_WW_%s_mj"%self.channel).getError()**2+
                                               self.workspace4fit_.var("rrv_number_WZ_%s_mj"%self.channel).getError()**2+
                                               self.workspace4fit_.var("rrv_number_WJets0_%s_mj"%self.channel).getError()**2));
        getattr(self.workspace4fit_,"import")(rrv_number_data_mj);
        
        ## make the final plot
        mplot = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(rrv_mass_j.getBins()));
        rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

        ## plot solid style 
        model_data.plotOn(mplot,RooFit.Name("WZ"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_WW_%s_mj,model_WZ_%s_mj"%(self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WZ"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
        model_data.plotOn(mplot,RooFit.Name("WW"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_WW_%s_mj"%(self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
        model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj"%(self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
        model_data.plotOn(mplot,RooFit.Name("STop"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj"%(self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
        model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model_WJets0_%s_mj"%self.channel),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));

            
        ### solid line
        model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj"%self.channel), RooFit.LineColor(kBlack), RooFit.LineWidth(2));
        model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj"%(self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2));
        model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj"%(self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2));
        model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_WW_%s_mj"%(self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(1));
        model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_WW_%s_mj,model_WZ_%s_mj"%(self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(1));


        rdataset_data_mj.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

        ### signal window zone with vertical lines
        lowerLine         = TLine(65,0.,65,mplot.GetMaximum()); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(7);
        upperLine         = TLine(105,0.,105,mplot.GetMaximum()); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(7);
        mplot.addObject(lowerLine);
        mplot.addObject(upperLine);
            
        ### legend of the plot
        self.leg = self.legend4Plot(mplot,0,1,0.3,0.,0.,0.,0,0,1);
        self.leg.SetHeader('pre-fit')
        mplot.addObject(self.leg);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.1);

        parameters_list = model_data.getParameters(rdataset_data_mj);
        
        mplot_pull = self.get_pull(rrv_mass_j,mplot);

        self.draw_canvas_with_pull( rrv_mass_j,rdataset_data_mj,mplot,mplot_pull,0,parameters_list,"%s/m_j_fitting/"%(self.plotsDir), "m_j_prefit",'',1,0,1)

    #################################################################################################
    #################################################################################################

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_mlvj_model_single_MC(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",in_file_name," ",label,"  ",mlvj_model,"  ",in_range," ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.channel+"_mlvj");
        rdataset.Print()
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,mlvj_model,"_mlvj",constrainslist,ismc);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.PrintLevel(-1) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.PrintLevel(-1), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
        self.fitresultsmlvj.append(rfresult)

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(rrv_mass_lvj.getBins()));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,6,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        nPar = rfresult.floatParsFinal().getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-nPar;
        #print mplot.chiSquare();
        #print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(nPar ,mplot.chiSquare(nPar)*ndof, ndof );
        datahist = rdataset.binnedClone( rdataset.GetName()+"_binnedClone",rdataset.GetName()+"_binnedClone" )
        rdataset.Print()
        
        ## get the pull 
        
        mplot_pull      = self.get_pull(rrv_mass_lvj,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
        
        self.draw_canvas_with_pull( rrv_mass_lvj, datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, logy);
        
        #@# make missing plot
        if 'TTBAR' in in_file_name and 'sb' in in_range:
                self.draw_canvas_with_pull( rrv_mass_lvj, datahist,mplot, mplot_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, 1)

         
        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_mlvj"%(label,in_range,self.channel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.channel+"_mlvj");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.channel+"_"+self.wtagger_label+"_mlvj_13TeV",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
            print "##################### diagonalize ";
            model_pdf_deco = Deco.diagonalize(model_pdf); ## diagonalize            
            print "##################### workspace for decorrelation ";
            wsfit_tmp.Print("v");
            print "##################### original  parameters ";
            model_pdf.getParameters(rdataset).Print("v");
            print "##################### original  decorrelated parameters ";
            model_pdf_deco.getParameters(rdataset).Print("v");
            print "##################### original  pdf ";
            model_pdf.Print();
            print "##################### decorrelated pdf ";
            model_pdf_deco.Print();

            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4fit_,"import")(model_pdf_deco);

            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame( RooFit.Bins(rrv_mass_lvj.getBins()));
            
            if label=="_TTbar" and in_range=="_sig":
                
                rdataset.plotOn(mplot_deco, RooFit.Name("Powheg Sample"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name("TTbar_Powheg"),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset = RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## draw the error band with the area
                self.workspace4fit_.var("rrv_number_TTbar_sig_%s_mlvj"%(self.channel)).Print();
            else:
                rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.leg = self.legend4Plot(mplot_deco,0); ## add the legend                
            mplot_deco.addObject(self.leg);

            self.draw_canvas( mplot_deco, "%s/other/"%(self.plotsDir), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print();

    #################################################################################################
    #################################################################################################

    ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

        print "############### Fit mj single MC sample",in_file_name," ",label,"  ",in_model_name," ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.channel+"_mj");
        rdataset_mj.Print();

        ## make the extended model
        model = self.make_Model(label,in_model_name);
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.Extended(kTRUE), RooFit.PrintLevel(-1) );
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.PrintLevel(-1), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
        getattr(self.workspace4fit_,'import')(rfresult)
        self.fitresultsmj.append(rfresult)
        
        ## Plot the result
        mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(rrv_mass_j.getBins()));
        rdataset_mj.plotOn( mplot, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,6,"L");
        ## re-draw the dataset
        rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## draw the function
        model.plotOn( mplot );# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = self.get_pull(rrv_mass_j, mplot);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        nPar = rfresult.floatParsFinal().getSize();
        #nBinX = mplot.GetNbinsX();
        nBinX = 22
        ndof  = nBinX-nPar;
        #print "#################### nPar=%s, nBinX=%s , chiSquare=%s/%s"%(nPar,nBinX,mplot.chiSquare('model%s_%s_mj_Norm[rrv_mass_j]'%(label,self.channel),'h_rdataset4fit%s_%s_mj'%(label,self.channel),nPar)*ndof,ndof);
        #print self.calculate_chi2(rdataset_mj,rrv_mass_j, mplot, ndof,1)
        mplot.Print()

        datahist = rdataset_mj.binnedClone( rdataset_mj.GetName()+"_binnedClone",rdataset_mj.GetName()+"_binnedClone" )
        parameters_list = model.getParameters(rdataset_mj);

        self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_j_fitting/"%(self.plotsDir), label+in_file_name, in_model_name)

        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print();
        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();


        fullInt           = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        signalInt         = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("sig"));
        fullInt_val        = fullInt.getVal()
        signalInt_val         = signalInt.getVal()/fullInt_val
        
        rrv_number_sig        = RooRealVar("rrv_number"+label+"_"+self.channel+"_sig","rrv_number"+label+"_"+self.channel+"_sig",self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal() * signalInt_val)
        getattr(self.workspace4limit_,'import')(rrv_number_sig)
        
        ##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV 
        par=parameters_list.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                #param.Print();
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
            param=par.Next()

    #################################################################################################
    #################################################################################################

    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label):# to get the shape of m_lvj,jet_mass="jet_mass_pr"

        if not options.read_trees:
                self.get_mj_and_mlvj_dataset_from_file(in_file_name, label)
        else:
                print "################### get_mj_and_mlvj_dataset : ",in_file_name,"  ",label,"  ##################";

                fileIn_name = TString(options.inPath+"/"+self.file_Directory+in_file_name);
                fileIn = TFile(fileIn_name.Data());
                treeIn = fileIn.Get("BasicTree");
                
                rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
                rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
                rrv_weight   = RooRealVar("rrv_weight","rrv_weight",-1000. ,10000000.)
                rrv_mass_j_sb_lo        = rrv_mass_j.clone("mj_sb_lo")
                rrv_mass_j_sb_hi        = rrv_mass_j.clone("mj_sb_hi")
                rrv_mass_j_sig                = rrv_mass_j.clone("mj_sig")
                rrv_mass_j_sb_lo.setRange(self.mj_sideband_lo_min,self.mj_sideband_lo_max)
                rrv_mass_j_sb_lo.setBins((self.mj_sideband_lo_max-self.mj_sideband_lo_min)/5)
                rrv_mass_j_sb_hi.setRange(self.mj_sideband_hi_min,self.mj_sideband_hi_max)
                rrv_mass_j_sb_hi.setBins((self.mj_sideband_hi_max-self.mj_sideband_hi_min)/5)
                rrv_mass_j_sig.setRange(self.mj_signal_min,self.mj_signal_max)
                rrv_mass_j_sig.setBins((self.mj_signal_max-self.mj_signal_min)/5)

                nbins_mlvj=int((rrv_mass_lvj.getMax()-rrv_mass_lvj.getMin())/self.BinWidth_mlvj);
                rrv_mass_lvj.setBins(nbins_mlvj);
                   
                ##### dataset of m_j -> scaled and not scaled to lumi 
                rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
                rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
                ##### dataset of m_lvj -> scaled and not scaled to lumi in different regions
                rdataset_sb_mlvj = RooDataSet("rdataset"+label+"_sb"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

                rdataset_sb_lo_mlvj = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

                rdataset_sig_mlvj = RooDataSet("rdataset"+label+"_sig"+"_"+self.channel+"_mlvj","rdataset"+label+"_sig"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

                rdataset_sb_hi_mlvj = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );
                

                rdataset4fit_sb_mlvj = RooDataSet("rdataset4fit"+label+"_sb"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );
                
                rdataset4fit_sb_lo_mlvj = RooDataSet("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );
                
                rdataset4fit_sig_mlvj = RooDataSet("rdataset4fit"+label+"_sig"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sig"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

                rdataset4fit_sb_hi_mlvj = RooDataSet("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );
                
                if True:#'data' in label:
                  ### datasets for simultaneous fit with combine
                  dataset_mj_sb_lo        = RooDataSet('dataset_2d_sb_lo_%s'%self.channel,'dataset_2d_sb_lo_%s'%self.channel,RooArgSet(rrv_mass_j_sb_lo, rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))
                  dataset_mj_sb_hi        = RooDataSet('dataset_2d_sb_hi_%s'%self.channel,'dataset_2d_sb_hi_%s'%self.channel,RooArgSet(rrv_mass_j_sb_hi, rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))
                  dataset_mj_sig        = RooDataSet('dataset_2d_sig_%s'%self.channel,'dataset_2d_sig_%s'%self.channel,RooArgSet(rrv_mass_j_sig, rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))


                ### categorize the event in sideband and signal region --> combined dataset 

                data_category = RooCategory("data_category","data_category");
                data_category.defineType("sideband");
                data_category.defineType("sig");
                combData = RooDataSet("combData"+label+"_"+self.channel,"combData"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
                combData4fit = RooDataSet("combData4fit"+label+"_"+self.channel,"combData4fit"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
                
                print "###### N entries: ", treeIn.GetEntries()


                for i in range(treeIn.GetEntries()):

                    if i % 10000 == 0: print "iEntry: ",i
                    treeIn.GetEntry(i);
                            
                    if i==0:
                        if TString(label).Contains('data'):
                            tmp_scale_to_lumi = 1
                        else:
                            #used to undo the scaling to lumi, needed for the fits to MC
                            tmp_scale_to_lumi = treeIn.totEventWeight / (treeIn.puweight*(treeIn.genweight/abs(treeIn.genweight)))

                    tmp_jet_mass=getattr(treeIn, 'Mjpruned');

                    self.isGoodEvent = 0;   
                    if treeIn.MWW> rrv_mass_lvj.getMin() and treeIn.MWW<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax():
                        self.isGoodEvent = 1;  
         
                    if self.isGoodEvent == 1:
                        ### weigh MC events                                                        
                        #@#change weight for new trees, lumi is already taken care of in totEventWeight
                        if TString(label).Contains('data'):
                            tmp_event_weight=1.;
                            tmp_event_weight4fit=1.; 
                        else:
                            tmp_event_weight                 = treeIn.totEventWeight
                            tmp_event_weight4fit         = tmp_event_weight/tmp_scale_to_lumi
                        #@#        
 
                        rrv_mass_lvj.setVal(treeIn.MWW);
                        #sideband lo only
                        if (tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max):
                            rdataset_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                            rdataset4fit_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                            
                            rdataset_sb_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                            rdataset4fit_sb_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

                            data_category.setLabel("sideband");
                            combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                            combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                                    
                        #signal region, blind data if needed
                        if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                            rdataset_sig_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                            rdataset4fit_sig_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                            
                            data_category.setLabel("sig");
                            combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                            combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

                        #sideband hi only, for sim-fit
                        if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                            rdataset_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                            rdataset4fit_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                            
                            rdataset_sb_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                            rdataset4fit_sb_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                            
                            data_category.setLabel("sideband");
                            combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                            combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);            
                            
                        rrv_mass_j.setVal( tmp_jet_mass );
                        #mj spectrum, cut out data in signal region if needed
                        if True:#'data' in label:
                          rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight )
                          rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit )
                          #datasets for simultaneaous fit in combine
                          ##sideband lo only
                          if tmp_jet_mass>= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
                            rrv_mass_j_sb_lo.setVal(tmp_jet_mass)
                            dataset_mj_sb_lo.add(RooArgSet(rrv_mass_j_sb_lo,rrv_mass_lvj), 1)
                          ##sideband hi
                          if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                            rrv_mass_j_sb_hi.setVal(tmp_jet_mass)
                            dataset_mj_sb_hi.add(RooArgSet(rrv_mass_j_sb_hi,rrv_mass_lvj), 1)        
                          ##signal
                          if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                            rrv_mass_j_sig.setVal(tmp_jet_mass)
                            dataset_mj_sig.add(RooArgSet(rrv_mass_j_sig,rrv_mass_lvj), 1)        


                ### scale to lumi for MC in 4fit datasets
                rrv_scale_to_lumi=RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi)
                rrv_scale_to_lumi.Print()
                getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)

                ### prepare m_lvj dataset to be compared with the fit results
                rrv_number_dataset_sig_mlvj=RooRealVar("rrv_number_dataset_sig"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_sig"+label+"_"+self.channel+"_mlvj",rdataset_sig_mlvj.sumEntries());
                rrv_number_dataset_AllRange_mlvj=RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj",rdataset_sig_mlvj.sumEntries()+rdataset_sb_mlvj.sumEntries());
                
                getattr(self.workspace4fit_,"import")(rrv_number_dataset_sig_mlvj)
                getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj)

                ### import the dataset       
                getattr(self.workspace4fit_,"import")(rdataset_sb_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset_sig_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset4fit_sb_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset4fit_sig_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj);
                getattr(self.workspace4fit_,"import")(rdataset_mj)
                getattr(self.workspace4fit_,"import")(rdataset4fit_mj)
                getattr(self.workspace4fit_,"import")(combData);
                getattr(self.workspace4fit_,"import")(combData4fit);
                if 'data' in label:
                    getattr(self.workspace4fit_,'import')(dataset_mj_sb_lo)
                    getattr(self.workspace4fit_,'import')(dataset_mj_sb_hi)
                    getattr(self.workspace4fit_,'import')(dataset_mj_sig)

                ### write in the output 
                self.file_out.write("\n%s events number in m_lvj from dataset: %s"%(label,rdataset_sig_mlvj.sumEntries()))

                #@#write datasets to file
                if 'WJets0_' in label:
                        fileOut        = TFile.Open('cards_%s_%s_900_3500/datasets_%s_%s.root'%(self.channel,self.wtagger_label,self.channel,self.wtagger_label),'recreate')
                else:
                        fileOut        = TFile.Open('cards_%s_%s_900_3500/datasets_%s_%s.root'%(self.channel,self.wtagger_label,self.channel,self.wtagger_label),'update')
                self.workspace4fit_.Write()
                fileOut.Close()
                #@#

                #### print everything
                
                rdataset_sb_lo_mlvj.Print();
                rdataset_sig_mlvj.Print();
                rdataset_sb_hi_mlvj.Print();
                rdataset_mj.Print();
                rdataset4fit_sb_lo_mlvj.Print();
                rdataset4fit_sig_mlvj.Print();
                rdataset4fit_sb_hi_mlvj.Print();
                rdataset4fit_mj.Print();
                combData.Print();
                combData4fit.Print();

    #################################################################################################
    #################################################################################################

    ###get MC and data histograms from file rather than reading TTrees
    def get_mj_and_mlvj_dataset_from_file(self,in_file_name, label):
        fileIn        = TFile.Open('cards_%s_%s_%s_%s/datasets_%s_%s.root'%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi,self.channel,self.wtagger_label))
        w_tmp        = fileIn.Get('workspace4fit_')
        getattr(self.workspace4fit_,'import')(w_tmp.var('rrv_scale_to_lumi'+label+'_'+self.channel))
        getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_sig"+label+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_sb"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_sig"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_sb"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_sig"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("combData"+label+"_"+self.channel))
        getattr(self.workspace4fit_,'import')(w_tmp.data("combData4fit"+label+"_"+self.channel))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_"+self.channel+"_mj"))
        getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_"+self.channel+"_mj"))
        if 'data' in label:
          getattr(self.workspace4fit_,'import')(w_tmp.data("dataset_2d_sb_lo_%s"%self.channel))
          getattr(self.workspace4fit_,'import')(w_tmp.data("dataset_2d_sb_hi_%s"%self.channel))
          getattr(self.workspace4fit_,'import')(w_tmp.data("dataset_2d_sig_%s"%self.channel))
        fileIn.Close()

    #################################################################################################
    #################################################################################################
        
    #### function to run the selection on data to build the datasets 
    def get_data(self):
        print "############### get_data ########################"
        self.get_mj_and_mlvj_dataset(self.file_data,"_data")
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_dataset_sig_data_%s_mlvj"%(self.channel)).clone("observation_for_counting"))

    #################################################################################################
    #################################################################################################

    #### Define the steps to fit STop MC in the mj and mlvj spectra
    def fit_STop(self):
        print "############################## fit_STop  #################################"
        self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")
        self.fit_mj_single_MC(self.file_STop_mc,"_STop","ExpGaus");
        self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop","_sb","Exp", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop","_sig","Exp", 0, 0, 1);
        print "________________________________________________________________________"

    #################################################################################################
    #################################################################################################

    ##### Define the steps to fit WW MC in the mj and mlvj spectra
    def fit_WW(self):
        print "############################# fit_WW ################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WW_mc,"_WW")
        self.fit_mj_single_MC(self.file_WW_mc,"_WW","2GausWW");       
        self.fit_mlvj_model_single_MC(self.file_WW_mc,"_WW","_sb","Exp", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_WW_mc,"_WW","_sig",self.MODEL_4_mlvj, 0, 0, 1);
        print "________________________________________________________________________"

    #################################################################################################
    #################################################################################################

    ##### Define the steps to fit WZ MC in the mj and mlvj spectra
    def fit_WZ(self):
        print "############################# fit_WZ ################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WZ_mc,"_WZ")
        self.fit_mj_single_MC(self.file_WZ_mc,"_WZ","2GausWZ");       
        self.fit_mlvj_model_single_MC(self.file_WZ_mc,"_WZ","_sb","ExpN", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_WZ_mc,"_WZ","_sig","ExpN", 0, 0, 1);
        print "________________________________________________________________________"   
        
    #################################################################################################
    #################################################################################################        
        
    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_TTbar(self):
        print "################################ fit_TTbar #########################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj
        self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar","2Gaus_ErfExp");
        self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar","_sb","ExpTail",1);
        self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar","_sig","ExpTail",1, 0, 1);
        print "________________________________________________________________________"

    #################################################################################################
    #################################################################################################

    ##### Define the steps to fit WJets MC in the mj and mlvj spectra
    def fit_WJets(self):
        print "######################### fit_WJets ########################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj
        #self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets1")# to get the shape of m_lvj
        ### Fit in mj depends on the mlvj lower limit -> fitting the turn on at low mass or not
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0","ErfExp")
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets1","User1");
        #### Fit the mlvj in sb_lo, signal region using two different model as done in the mj

        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sig",self.MODEL_4_mlvj, 0, 0, 1, 1);
        #self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets1","_sb",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);
        #self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets1","_sig",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);       
        print "________________________________________________________________________"

    #################################################################################################
    #################################################################################################

    ##### Analysis with sideband alpha correction 
    def analysis_sideband_correction_method1(self):

        print "##################### Start sideband correction full analysis ##############"

        ### Fit all MC components in both mj and mlvj
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
        self.fit_WJets()
        self.fit_TTbar()
        self.fit_WW()
        self.fit_WZ()
        self.fit_STop()
        print "________________________________________________________________________"

        ### take the real data
        self.get_data()
        ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
        self.fit_WJetsNorm();
        ### fit data in the mlvj low sideband with two different models
        #self.fit_mlvj_in_Mj_sideband("_WJets1","_sb",self.MODEL_4_mlvj_alter,1)
        self.fit_mlvj_in_Mj_sideband("_WJets0","_sb",self.MODEL_4_mlvj,1)

        ### Prepare the workspace and datacards     
        self.prepare_limit("sideband_correction_method1",1,0,0)
        ### finale plot and check of the workspace
        self.read_workspace(1)
        ### print all fitresults
        for results in [self.fitresultsmj,self.fitresultsmlvj,self.fitresultsfinal]:
            for i in results:
                i.Print()
        
    #################################################################################################
    #################################################################################################

    ##### Prepare the workspace for the limit and to store info to be printed in the datacard
    ###import and rename everything needed in the final limit setting procedure
    def prepare_limit(self,mode, isTTbarFloating=0, isVVFloating=0, isSTopFloating=0):

        if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail":
            self.workspace4fit_.var("Deco_WJets0_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
            self.workspace4fit_.var("Deco_WJets0_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
            self.workspace4fit_.var("Deco_WJets0_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
            self.workspace4fit_.var("Deco_WJets0_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);

            ### TTbar use expN
            if isTTbarFloating !=0:
                self.workspace4fit_.var("Deco_TTbar_sig_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                self.workspace4fit_.var("Deco_TTbar_sig_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                self.workspace4fit_.var("Deco_TTbar_sb_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                self.workspace4fit_.var("Deco_TTbar_sb_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);

        #@#prepare functions for simultaneous fit in combine
        #all pdfs are splitted into sb and sig regions (only in mj spectrum, mlvj is done in the signal code)
        getattr(self.workspace4fit_,'import')(self.workspace4fit_.pdf("model_pdf_WJets0_%s_mj"%self.channel).clone("model_pdf_WJets_%s_mj"%self.channel))
        bkgs                = ['WJets','TTbar','WW','WZ','STop']
        ras_mass_j        = RooArgSet(self.workspace4fit_.var('rrv_mass_j'))
        for bkg in bkgs:
            getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf('model_pdf_%s_%s_mj'%(bkg,self.channel)).clone('mj_%s_%s'%(bkg,self.channel)))
        for region in ['sig','sb_lo','sb_hi']:
            for bkg in bkgs:
                custom_mj        = RooCustomizer(self.workspace4fit_.pdf('model_pdf_%s_%s_mj'%(bkg,self.channel)),'%s_mj_%s'%(bkg,region))
                int_all                = self.workspace4fit_.pdf('model_pdf_%s_%s_mj'%(bkg,self.channel)).createIntegral(ras_mass_j,ras_mass_j,region).getVal()
                #need to rescale fractions (rrv_frac_...) for different regions
                if bkg=='TTbar':
                    old_frac        = self.workspace4fit_.var('rrv_frac_TTbar_'+self.channel)
                    new_frac        = RooRealVar('rrv_frac_'+region+'_'+bkg+'_'+self.channel,'rrv_frac_'+region+'_'+bkg+'_'+self.channel,0)
                    int_gaus        = self.workspace4fit_.pdf('gaus1_TTbar_'+self.channel).createIntegral(ras_mass_j,ras_mass_j,region).getVal()*(1-old_frac.getVal())
                    int_erfexp        = self.workspace4fit_.pdf('erfexp_TTbar_'+self.channel+'_mj').createIntegral(ras_mass_j,ras_mass_j,region).getVal()*old_frac.getVal()
                    new_frac.setVal(int_erfexp/(int_gaus+int_erfexp))
                    custom_mj.replaceArg(old_frac,new_frac)
                elif bkg=='STop':
                    old_frac        = self.workspace4fit_.var('rrv_high_STop_'+self.channel)
                    new_frac        = RooRealVar('rrv_frac_'+region+'_'+bkg+'_'+self.channel,'rrv_frac_'+region+'_'+bkg+'_'+self.channel,0)
                    int_gaus        = self.workspace4fit_.pdf('gaus_STop_'+self.channel).createIntegral(ras_mass_j,ras_mass_j,region).getVal()*(1-old_frac.getVal())
                    int_exp        = self.workspace4fit_.pdf('exp_STop_'+self.channel).createIntegral(ras_mass_j,ras_mass_j,region).getVal()*old_frac.getVal()
                    new_frac.setVal(int_exp/(int_gaus+int_exp))
                    custom_mj.replaceArg(old_frac,new_frac)
                elif bkg=='WW' or bkg=='WZ':
                    old_frac1        = self.workspace4fit_.var('rrv_frac1_%s_'%bkg+self.channel)
                    new_frac11        = RooRealVar('rrv_frac11_'+region+'_'+bkg+'_'+self.channel,'rrv_frac11_'+region+'_'+bkg+'_'+self.channel,0)
                    int_gaus1        = self.workspace4fit_.pdf('gaus1_%s_'%bkg+self.channel).createIntegral(ras_mass_j,ras_mass_j,region).getVal()*old_frac1.getVal()
                    int_gaus2        = self.workspace4fit_.pdf('gaus2_%s_'%bkg+self.channel).createIntegral(ras_mass_j,ras_mass_j,region).getVal()*(1-old_frac1.getVal())
                    new_frac11.setVal(((int_gaus1)/(int_gaus1+int_gaus2)))
                    custom_mj.replaceArg(old_frac1,new_frac11)
                custom_mj.replaceArg(self.workspace4fit_.var('rrv_mass_j'),self.workspace4fit_.var('mj_%s'%region))
                m_pruned_pdf        = custom_mj.build()
                m_pruned_pdf.Print()
                getattr(self.workspace4limit_,'import')(m_pruned_pdf.clone('%s_mj_%s_%s'%(bkg,region,self.channel)),RooFit.RecycleConflictNodes())
        
        self.workspace4fit_.allPdfs().Print("V")
        print "####################### prepare_limit for %s method ####################"%(mode);
        
        rrv_x = self.workspace4limit_.var("rrv_mass_lvj");

        #take decorellated ttbar
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_pdf_TTbar_sb_%s_mlvj_Deco_TTbar_sb_%s_HPV_mlvj_13TeV"%(self.channel,self.channel)).clone('TTbar_mlvj_sb_%s'%self.channel))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_pdf_TTbar_sig_%s_mlvj_Deco_TTbar_sig_%s_HPV_mlvj_13TeV"%(self.channel,self.channel)).clone('TTbar_mlvj_sig_%s'%self.channel))
        #normal pdfs for rest
        for label in ['STop','WW','WZ']:
            getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_pdf_"+label+"_sb_"+self.channel+'_mlvj').clone('%s_mlvj_sb_%s'%(label,self.channel)))
            getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_pdf_"+label+"_sig_"+self.channel+'_mlvj').clone('%s_mlvj_sig_%s'%(label,self.channel)))
        for label in ['TTbar','STop','WW','WZ']:
            self.fix_Pdf(self.workspace4limit_.pdf('%s_mlvj_sig_%s'%(label,self.channel)), RooArgSet(rrv_x) ); 
            getattr(self.workspace4limit_,'import')(self.workspace4fit_.var("rrv_number_"+label+"_"+self.channel+"_mj").clone('norm_%s_%s'%(label, self.channel)))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.var("rrv_number_WJets0_"+self.channel+"_mj").clone('norm_WJets_%s'%self.channel))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_pdf_WJets0_sb_"+self.channel+"_mlvj").clone("WJets_mlvj_sb_"+self.channel))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_pdf_WJets0_sig_%s_undeco_mlvj"%self.channel).clone("WJets_mlvj_sig_%s"%self.channel), RooFit.RecycleConflictNodes())
        self.workspace4limit_.pdf("WJets_mlvj_sig_%s"%self.channel).Print()
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf('model_WJets0_%s_mj'%self.channel).clone('WJets_mj_%s'%self.channel))

        #get m_pruned and mlvj datasets
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('rdataset_data_sb_lo_%s_mlvj'%self.channel))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('rdataset_data_sig_%s_mlvj'%(self.channel)).Clone('data_obs_%s_%s'%(self.channel,self.wtagger_label)))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('dataset_2d_sb_lo_%s'%self.channel))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('dataset_2d_sb_hi_%s'%self.channel))
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('dataset_2d_sig_%s'%self.channel))        
        getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('rdataset_data_%s_mj'%self.channel))

        if isTTbarFloating:
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_sig_%s_mlvj_Deco_TTbar_sig_%s_%s_mlvj_13TeV"%(self.channel, self.channel, self.wtagger_label)).clone("TTbar_%s_%s"%(self.channel,self.wtagger_label)))
        else :
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_sig_%s_mlvj"%(self.channel)).clone("TTbar_%s_%s"%(self.channel,self.wtagger_label)))

        print " ############## Workspace for limit ";
        parameters_workspace = self.workspace4limit_.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next()

        ### main modality for the alpha function method
        self.workspace4limit_.allVars().Print("V")
        self.workspace4fit_.allVars().Print("V")
        self.workspace4fit_.allPdfs().Print("V")

        ### Save the workspace
        self.save_workspace_to_file();

    #################################################################################################
    #################################################################################################

    #### Method used in order to save the workspace in a output root file
    def save_workspace_to_file(self):
        self.workspace4limit_.writeToFile(self.file_rlt_root);
        self.file_out.close()

    #################################################################################################
    #################################################################################################

    #### Read the final workspace and produce the latest plots
    def read_workspace(self, logy=0):

        ### Taket the workspace for limits  
        file = TFile(self.file_rlt_root) ;
        workspace = file.Get("workspace4limit_") ;

        workspace.data("data_obs_%s_%s"%(self.channel,self.wtagger_label)).Print()

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_%s_%s"%(self.channel,self.wtagger_label));

        model_pdf_WJets  = workspace.pdf("WJets_mlvj_sig_%s"%self.channel).clone("WJets");
        model_pdf_WW     = workspace.pdf("WW_mlvj_sig_%s"%self.channel).clone("WW");
        model_pdf_WZ     = workspace.pdf("WZ_mlvj_sig_%s"%self.channel).clone("WZ");
        model_pdf_TTbar  = workspace.pdf("TTbar_mlvj_sig_%s"%self.channel).clone("TTbar");
        model_pdf_STop   = workspace.pdf("STop_mlvj_sig_%s"%self.channel).clone("STop");
        
        model_pdf_WJets.Print();
        model_pdf_WW.Print();
        model_pdf_WZ.Print();
        model_pdf_TTbar.Print();
        model_pdf_STop.Print();
        
        rrv_number_WJets  = workspace.var("rrv_number_WJets0_"+self.channel+"_sig");
        rrv_number_WW     = workspace.var("rrv_number_WW_"+self.channel+"_sig");
        rrv_number_WZ     = workspace.var("rrv_number_WZ_"+self.channel+"_sig");
        rrv_number_TTbar  = workspace.var("rrv_number_TTbar_"+self.channel+"_sig");
        
        rrv_number_STop   = workspace.var("rrv_number_STop_"+self.channel+"_sig");

        rrv_number_WJets.Print();
        rrv_number_WW.Print();
        rrv_number_WZ.Print();
        rrv_number_TTbar.Print();
        rrv_number_STop.Print();
        
        #### Prepare the final plot starting from total background 
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC","rrv_number_Total_background_MC",
                rrv_number_WJets.getVal()+rrv_number_WW.getVal()+rrv_number_WZ.getVal()+rrv_number_TTbar.getVal()+rrv_number_STop.getVal());

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets.getError()**2+rrv_number_WW.getError()**2+rrv_number_WZ.getError()**2+rrv_number_TTbar.getError()**2+rrv_number_STop.getError()**2
                ));

        #### Total pdf 
        model_Total_background_MC = RooAddPdf("model_Total_background_MC","model_Total_background_MC",RooArgList(model_pdf_WJets,model_pdf_WW,model_pdf_WZ,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_WW,rrv_number_WZ,rrv_number_TTbar,rrv_number_STop));


        if data_obs.sumEntries() != 0:
         #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
         scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
        else:
         scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()   
        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(rrv_x.getBins()));
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        #plot pdfs
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets,WW,WZ,TTbar,STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar"), RooFit.Components("WW,WZ,TTbar,STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WW"), RooFit.Components("WW,WZ,STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW"]), RooFit.LineColor(kBlack), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WZ"), RooFit.Components("WZ,STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WZ"]), RooFit.LineColor(kBlack), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop"), RooFit.Components("STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());

        #solid line
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets,WW,WZ,TTbar,STop"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar_line_invisible"), RooFit.Components("WW,WZ,TTbar,STop"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WW_line_invisible"), RooFit.Components("WW,WZ,STop"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WZ_line_invisible"), RooFit.Components("WZ,STop"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop_line_invisible"), RooFit.Components("STop"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());


        #### plot the observed data using poissonian error bar
        self.getData_PoissonInterval(data_obs,mplot);
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Invisible());
        
        mplot_pull = self.get_pull(rrv_x,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC,self.FloatingParams,workspace ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        leg =        self.legend4Plot(mplot,0,1,0.25,0.,0.1,0.,0);
        leg.SetHeader('pre-fit')
        mplot.addObject(leg);
        
        mplot.GetYaxis().SetRangeUser(1e-1,5e3);
        mplot.GetXaxis().SetTitle('M_{WV} (GeV)')
            
        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal());
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-self.nPar_float_in_fitTo;
        #print "nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof );
        datahist = data_obs.binnedClone( data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone" )
        
        parameters_list = RooArgList();
        self.draw_canvas_with_pull( rrv_x,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir),"check_workspace_for_limit","",0,1,0,1);
    
    #################################################################################################
    #################################################################################################

### funtion to run the complete alpha analysis
def pre_limit_sb_correction(method, channel, in_mj_min=40, in_mj_max=150, in_mlvj_min=900, in_mlvj_max=3500, fit_model="ExpN", fit_model_alter="ExpTail"): 

    print "#################### pre_limit_sb_correction: channel %s, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    
    boostedW_fitter=doFit_wj_and_wlvj(channel,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    getattr(boostedW_fitter,"analysis_sideband_correction_%s"%(method) )();

    #################################################################################################
    #################################################################################################

#### Main Code
if __name__ == '__main__':

    channel=options.channel;
            
    pre_limit_sb_correction("method1",channel,40,150,options.mlvj_lo,3500,"ExpN","ExpTail")

