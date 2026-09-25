////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(6) spin         //
//                     in momentum space up to 21 orbitals                    //
//                with a cap on the number of particles per band              //
//                                                                            //
//                        last modification : 04/11/2024                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU6SPINANDCAPMOMENTUMSPACELONG_H
#define FERMIONONSQUARELATTICEWITHSU6SPINANDCAPMOMENTUMSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong.h"

#include <iostream>



class FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong : public FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong
{

 protected:

  // minimum number of particles in band 0
  int MinNbrParticlesBand0;
  // minimum number of particles in band 1
  int MinNbrParticlesBand1;
  // minimum number of particles in band 2
  int MinNbrParticlesBand2;

  // maximum number of particles in band 0
  int MaxNbrParticlesBand0;
  // maximum number of particles in band 1
  int MaxNbrParticlesBand1;
  // maximum number of particles in band 2
  int MaxNbrParticlesBand2;

 public:

  // default constructor
  //
  FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong();
  
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
  FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

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
  FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // constructor
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
  // totalSz = twice the total Sz (or any U(1) quantum number)
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, int totalSz, unsigned long memory);

  // copy constructor (without duplicating data)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong(const FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong ();

  // assignement (without duplicating data)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong& operator = (const FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpaceLong& fermions);

  // clone Hilbert space (without duplicating data)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension for a single band
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // singleBandTotalKx = total momentum along x
  // singleBandTotalKy = total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy);

  // evaluate Hilbert space dimension for a single band
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalSz = current total spin 
  // singleBandTotalKx = total momentum along x
  // singleBandTotalKy = total momentum along y
  // singleBandTotalSz = total spin 
  // return value = Hilbert space dimension
  virtual long EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentTotalSz, int singleBandTotalKx, int singleBandTotalKy, int singleBandTotalSz);

  // generate all states corresponding to the constraints for a single band
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // singleBandTotalKx = total momentum along x
  // singleBandTotalKy = total momentum along y
  // singleBandStateDescription = pointer to the single band state description array
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, ULONGLONG* singleBandStateDescription, long pos);
  
  // generate all states corresponding to the constraints for a single band
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalSz = current total spin 
  // singleBandTotalKx = total momentum along x
  // singleBandTotalKy = total momentum along y
  // singleBandTotalSz = total spin 
  // singleBandStateDescription = pointer to the single band state description array
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentTotalSz, int singleBandTotalKx, int singleBandTotalKy, int singleBandTotalSz, ULONGLONG* singleBandStateDescription, long pos);

  // evaluate all the single band Hilbert spaces
  //
  // maxBandOccupation = maiximum occupation of a single band
  // singleBandTotalKxMax = reference on the array for the maximum total Kx values 
  // singleBandTotalKyMax = reference on the array for the maximum total Ky values
  // singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky)
  // singleBandStates = reference on the array for all the Hilbert space basis states 
  virtual void GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, ULONGLONG****& singleBandStates);

  // evaluate all the single band Hilbert spaces
  //
  // maxBandOccupation = maximum occupation of a single band
  // singleBandTotalKxMax = reference on the array for the maximum total Kx values 
  // singleBandTotalKyMax = reference on the array for the maximum total Ky values
  // singleBandTotalSz = reference on the array for the maximum total Sz values
  // singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky, fourth index being the total Sz)
  // singleBandStates = reference on the array for all the Hilbert space basis states 
  virtual void GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, int*& singleBandTotalSz, long****& singleBandHilbertDimensions, ULONGLONG*****& singleBandStates);
  
  // generate all states using the single band hilbert spaces
  //
  virtual void GenerateStatesFromSingleBandHilbertSpaces();
  
  // generate all states using the single band hilbert spaces for the spin conserved case
  //
  virtual void GenerateSpinConvervedStatesFromSingleBandHilbertSpaces();


};


#endif


