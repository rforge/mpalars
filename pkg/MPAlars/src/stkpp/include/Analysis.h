/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project:  stkpp::Analysis
 * Purpose:  Primary include file for Analysis project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file Analysis.h
 *  @brief This file include all the header files of the project Analysis.
 * 
 *  @defgroup Analysis Special functions tools
 *  @brief In this project we compute usual and special functions.
 *  The Analysis provide al the tools necessary to the computation of the usual
 *  and special functions like gamma, beta, gamma ratio, beta ration functions.
 *  It provide generic algorithms and the usual mathematical constants.
 **/

/**
 * @ingroup Analysis
 * @namespace STK::Const
 * @brief This is the namespace enclosing the usual mathematical constants
 **/
/**
 * @ingroup Analysis
 * @namespace STK::Funct
 * @brief The namespace Funct enclose all usual and special functions.
 * The namespace Funct is the domain space of the special function
 * like gamma function, beta function, incomplete gamma function,
 * incomplete beta function... It include also some useful raw
 * functions like log1p...
 **/

#ifndef ANALYSIS_H
#define ANALYSIS_H

// templated generic algorithms
#include "../projects/Analysis/include/STK_Algo.h"

// namespace Const
// Mathematical constant
#include "../projects/Analysis/include/STK_Const_Math.h"

// namespace Funct
// usual fonctions
#include "../projects/Analysis/include/STK_Funct_util.h"
// raw functions
#include "../projects/Analysis/include/STK_Funct_raw.h"
// gamma function
#include "../projects/Analysis/include/STK_Funct_gamma.h"
// gamma ratio function
#include "../projects/Analysis/include/STK_Funct_gammaRatio.h"
// beta Ratio function
#include "../projects/Analysis/include/STK_Funct_betaRatio.h"


#endif /*ANALYSIS_H*/
