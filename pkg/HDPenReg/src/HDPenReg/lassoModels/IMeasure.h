/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 5 nov. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file IMeasure.h
 *  @brief In this file .
 **/


#ifndef IMEASURE_H_
#define IMEASURE_H_

#include "../../stkpp/include/STKpp.h"

class IMeasure
{
  public:
    IMeasure(){};
    virtual ~IMeasure(){};

    virtual STK::Real measure(STK::CVectorX const& yTrue, STK::CVectorX &yEstimate) {return 0.;};
};


class Residuals : public IMeasure
{
  public:
    Residuals() : IMeasure(){};

    STK::Real measure(STK::CVectorX const& yTrue, STK::CVectorX &yEstimate)
    {
      STK::Real res;
      res=(yTrue-yEstimate).norm2()/(yTrue.sizeRows());

      return res;
    }
};

class AUC : public IMeasure
{
  public:
    AUC() : IMeasure(){};

    STK::Real measure(STK::CVectorX const& yTrue, STK::CVectorX &yEstimate)
    {
      STK::Real auc, nb1, nb0, rsum1 = 0;
      nb1 = yTrue.sum();
      nb0 = yTrue.sizeRows() - nb1;

      STK::Array2DVector<int> index(yTrue.sizeRows());
      for(int i = 1; i <= index.sizeRows(); i++)
        index[i] = i;

      STK::Array2DVector<int> index2(index);

      int i, j,elem2;
      STK::Real elem1;

      for (i = 2; i <= yEstimate.sizeRows(); ++i)
      {
          elem1 = yEstimate[i];
          elem2 = index[i];
          for (j = i; j > 1 && yEstimate[j-1] > elem1; j--)
          {
            yEstimate[j] = yEstimate[j-1];
            index[j] = index[j-1];
          }
          yEstimate[j] = elem1;
          index[j] = elem2;
      }

      for (i = 2; i <= index.sizeRows(); ++i)
      {
          elem1 = index[i];
          elem2 = index2[i];
          for (j = i; j > 1 && index[j-1] > elem1; j--)
          {
            index[j] = index[j-1];
            index2[j] = index2[j-1];
          }
          index[j] = elem1;
          index2[j] = elem2;
      }
      for(i = 1; i <= yTrue.sizeRows(); i++)
      {
        if(yTrue[i]==1)
          rsum1 += index2[i];
      }

      auc = (rsum1-nb1*(nb1+1.)/2.)/(nb1 * nb0);

      return auc;
    }
};

#endif /* IMEASURE_H_ */
