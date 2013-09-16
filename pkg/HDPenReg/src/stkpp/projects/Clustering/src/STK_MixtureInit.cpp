/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * Project:  stkpp::Clustering
 * created on: 24 août 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureInit.cpp
 *  @brief In this file we implement the initialization methods.
 **/

#include "../include/STK_MixtureInit.h"
#include "../include/STK_IMixtureModelBase.h"

namespace STK
{

/* call the classInit() model initialization.
 * @return @c true if no error occur, @c false otherwise*/
bool RandomInit::run()
{
  try
  {
    p_model_->randomInit();
  }
  catch (Exception const& e)
  {
     msg_error_ = e.error();
     return false;
  }
  return true;
}

/* call the classInit() model initialization.
 * @return @c true if no error occur, @c false otherwise*/
bool ClassInit::run()
{
  try
  {
    p_model_->classInit();
  }
  catch (Exception const& e)
  {
     msg_error_ = e.error();
     return false;
  }
  return true;
}

/* call the classInit() model initialization.
 * @return @c true if no error occur, @c false otherwise*/
bool FuzziInit::run()
{
  try
  {
    p_model_->fuzziInit();
  }
  catch (Exception const& e)
  {
     msg_error_ = e.error();
     return false;
  }
  return true;
}

} // namespace STK
