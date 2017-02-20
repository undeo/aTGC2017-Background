/* specfunc/gsl_sf_hyperg.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_HYPERG_H__
#define __GSL_SF_HYPERG_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

namespace Mathelp {


/* Confluent hypergeometric function.
 * 1F1[a,b,x] = M(a,b,x)
 *
 * exceptions:
 */
int gsl_sf_hyperg_1F1_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_hyperg_1F1(double a, double b, double x);

/* Gauss hypergeometric function 2F1[a,b,c,x]
 * |x| < 1
 *
 * exceptions:
 */
int mgsl_sf_hyperg_2F1_e(double a, double b, const double c, const double x, gsl_sf_result * result);
double mgsl_sf_hyperg_2F1(double a, double b, double c, double x);


/* Gauss hypergeometric function
 * 2F1[aR + I aI, aR - I aI, c, x]
 * |x| < 1
 *
 * exceptions:
 */
int mgsl_sf_hyperg_2F1_conj_e(const double aR, const double aI, const double c, const double x, gsl_sf_result * result);
double mgsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);


/* Renormalized Gauss hypergeometric function
 * 2F1[a,b,c,x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int mgsl_sf_hyperg_2F1_renorm_e(const double a, const double b, const double c, const double x, gsl_sf_result * result);
double mgsl_sf_hyperg_2F1_renorm(double a, double b, double c, double x);


/* Renormalized Gauss hypergeometric function
 * 2F1[aR + I aI, aR - I aI, c, x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int mgsl_sf_hyperg_2F1_conj_renorm_e(const double aR, const double aI, const double c, const double x, gsl_sf_result * result);
double mgsl_sf_hyperg_2F1_conj_renorm(double aR, double aI, double c, double x);


}
__END_DECLS

#endif /* __GSL_SF_HYPERG_H__ */
