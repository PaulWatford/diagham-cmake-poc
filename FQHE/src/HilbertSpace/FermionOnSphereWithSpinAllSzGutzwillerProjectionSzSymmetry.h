/////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                       all Sz sectors, Sz<->-Sz symmetry                    //
//                  and using Gutziller projection in orbital space           //
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


#ifndef FERMIONONSPHEREWITHSPINALLSZGUTZWILLERPROJECTIONSZSYMMETRY_H
#define FERMIONONSPHEREWITHSPINALLSZGUTZWILLERPROJECTIONSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



class FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry :  public FermionOnSphereWithSpinAllSzSzSymmetry
{


  friend class FermionOnSphereWithSpinAllSzLzSzSymmetry;
  
 protected:


 public:

  // default constructor 
  //
  FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // minusParity = select the Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry (int nbrFermions, int totalLz, int lzMax, bool minusParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry(const FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry& operator = (const FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry& fermions);

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


