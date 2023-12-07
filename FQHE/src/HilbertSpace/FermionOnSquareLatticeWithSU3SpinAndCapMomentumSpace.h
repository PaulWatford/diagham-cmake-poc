////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(3) spin         //
//       in momentum space with a cap on the number of particles per band     //
//                                                                            //
//                        last modification : 03/12/2023                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU3SPINANDCAPMOMENTUMSPACE_H
#define FERMIONONSQUARELATTICEWITHSU3SPINANDCAPMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"

#include <iostream>



class FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace : public FermionOnSquareLatticeWithSU3SpinMomentumSpace
{

 protected:

  // maximum number of particles in band 0
  int MaxNbrParticlesBand0;
  // maximum number of particles in band 1
  int MaxNbrParticlesBand1;
  // maximum number of particles in band 2
  int MaxNbrParticlesBand2;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // maxNbrParticlesBand0 = maximum number of particles in band 0
  // maxNbrParticlesBand1 = maximum number of particles in band 1
  // maxNbrParticlesBand2 = maximum number of particles in band 2
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace(const FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace& operator = (const FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();


 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // maxNbrParticlesBand0 = current maximum number of particles in band 0
  // maxNbrParticlesBand1 = current maximum number of particles in band 1
  // maxNbrParticlesBand2 = current maximum number of particles in band 2
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // maxNbrParticlesBand0 = current maximum number of particles in band 0
  // maxNbrParticlesBand1 = current maximum number of particles in band 1
  // maxNbrParticlesBand2 = current maximum number of particles in band 2
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, long pos);


};


#endif


