////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(3) spin         //
//                  for more than 21 orbitals in momentum space               //
//                       with an upper and lower bounds                       //
//                     on the number of particles per band                    //
//                                                                            //
//                        last modification : 18/04/2024                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU3SPINANDMINMAXCAPMOMENTUMSPACELONG_H
#define FERMIONONSQUARELATTICEWITHSU3SPINANDMINMAXCAPMOMENTUMSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong.h"

#include <iostream>



class FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong : public FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong
{

 protected:

  // maximum number of particles in band 0
  int MinNbrParticlesBand0;
  // maximum number of particles in band 1
  int MinNbrParticlesBand1;
  // maximum number of particles in band 2
  int MinNbrParticlesBand2;

 public:

  // default constructor
  //
  FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong();
  
  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // minNbrParticlesBand0 = minimum number of particles in band 0
  // minNbrParticlesBand1 = minimum number of particles in band 1
  // minNbrParticlesBand2 = minimum number of particles in band 2
  // maxNbrParticlesBand0 = maximum number of particles in band 0
  // maxNbrParticlesBand1 = maximum number of particles in band 1
  // maxNbrParticlesBand2 = maximum number of particles in band 2
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong(const FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong& operator = (const FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong& fermions);

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


