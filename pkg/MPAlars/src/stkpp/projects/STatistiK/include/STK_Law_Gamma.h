/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2008  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::STatistiK::Law
 * Purpose:  Cauchy probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_LawGamma.h
 *  @brief In this file we define the Gamma probability distribution.
 **/

#ifndef STK_LAW_GAMMA_H
#define STK_LAW_GAMMA_H

#include "STK_Law_IUnivLaw.h"

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Gamma distribution law.
 * 
 *  The Gamma distribution is a continuous probability distribution
 *  The probability pdf function of the gamma distribution
 *  can be expressed in terms of the @ref gamma function:
 *  \f[
 *   f(x;k,\theta) = \left(\frac{x}{\theta}\right)^{\alpha-1}
 *                   \frac{e^{-x/\theta}}{\theta \, \Gamma(\alpha)} 
 *   \ \mathrm{for}\ x > 0 
 *  \f]
 *  where @em k>0 is the shape parameter
 *  and   \f$ \theta > 0 \f$ is the scale parameter.
 **/
class Gamma : public IUnivLaw<Real>
{
  public:
    /** Default constructor.
     *  @param shape shape (position) parameter
     *  @param scale scale (dispersion) parameter
     **/
    Gamma( const Real& shape, const Real& scale);
    /** Destructor. */
	  virtual ~Gamma();
    /** @return a pseudo Gamma random variate. */
    virtual Real rand() const;
    /** @return the value of the pdf
     *  @param x a positive real value
     **/
    virtual Real pdf(const Real& x) const;
    /** @return the value of the log-pdf
     *  @param x a positive real value
     **/
    virtual Real lpdf(const Real& x) const;
    /** @return the cumulative distribution function
     *  @param t a positive real value
     **/
    virtual Real cdf(const Real& t) const;
    /** @return the inverse cumulative distribution function
     *  @param p a probability number
     **/
    virtual Real icdf(const Real& p) const;
    /** @return a pseudo Gamma random variate with the specified parameters.
     *  @param shape shape (position) parameter
     *  @param scale scale (dispersion) parameter
     **/
    static Real rand( const Real& shape, const Real& scale);
    /** @return the value of the pdf
     *  @param x a positive real value
     **/
    static Real pdf(const Real& x, const Real& shape, const Real& scale);
    /** @return the value of the log-pdf
     *  @param x a positive real value
     **/
    static Real lpdf(const Real& x, const Real& shape, const Real& scale);

  protected:
    /** The shape parameter */
    Real a_;
    /** The scale parameter */
    Real b_;
    /** First and second constants for rand */
    mutable Real c_, d_;

};

} // namespace Law

} // namespace STK

#endif /*STK_LAWGAMMA_H*/
