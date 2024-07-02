////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                     without Sz conservation and using Gutziller            //
//                          projection in orbital space                       //
//                                                                            //
//                        last modification : 08/12/2023                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef FERMIONONSPHEREWITHSPINALLSZGUTZWILLERPROJECTION_H
#define FERMIONONSPHEREWITHSPINALLSZGUTZWILLERPROJECTION_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereWithSpinAllSzGutzwillerProjection :  public FermionOnSphereWithSpinAllSz
{

  friend class FermionOnSphereWithSpinAllSzLzSymmetry;
  friend class FermionOnSphereWithSpinAllSzLzSzSymmetry; 
  friend class FermionOnSphereWithSpinAllSzSzSymmetry; 
/*   friend class FermionOnSphereWithSpinAllSzLzSymmetry; */

 protected:


 public:

  // default constructor
  //
  FermionOnSphereWithSpinAllSzGutzwillerProjection();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinAllSzGutzwillerProjection (int nbrFermions, int totalLz, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinAllSzGutzwillerProjection(const FermionOnSphereWithSpinAllSzGutzwillerProjection& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinAllSzGutzwillerProjection ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinAllSzGutzwillerProjection& operator = (const FermionOnSphereWithSpinAllSzGutzwillerProjection& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a given state from a generic basis to the gutzwiller basis
  //
  // state = reference on the vector to convert
  // basis = pointer to the basis associated to state
  // return value = converted vector
  virtual RealVector GutzwillerProjection(RealVector& state, ParticleOnSphereWithSpin* basis);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);


  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up ( omitted: int totalSpin)
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, long pos);
  
};

#endif


