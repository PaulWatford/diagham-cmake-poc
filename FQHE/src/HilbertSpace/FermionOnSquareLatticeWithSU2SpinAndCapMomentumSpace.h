////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(2) spin         //
//       in momentum space with a cap on the number of particles per band     //
//                                                                            //
//                        last modification : 18/12/2024                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU2SPINANDCAPMOMENTUMSPACE_H
#define FERMIONONSQUARELATTICEWITHSU2SPINANDCAPMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include <iostream>



class FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace : public FermionOnSquareLatticeWithSU2SpinMomentumSpace
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:

  // minimum number of particles in band 0
  int MinNbrParticlesBand0;
  // minimum number of particles in band 1
  int MinNbrParticlesBand1;

  // maximum number of particles in band 0
  int MaxNbrParticlesBand0;
  // maximum number of particles in band 1
  int MaxNbrParticlesBand1;

public:

  // default constructor
  //
  FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace();
  
  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // maxNbrParticlesBand0 = maximum number of particles in band 0
  // maxNbrParticlesBand1 = maximum number of particles in band 1
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int kxMomentum, int kyMomentum, char* outputDirectory = 0, unsigned long memory = 10000000);

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // minNbrParticlesBand0 = minimum number of particles in band 0
  // minNbrParticlesBand1 = minimum number of particles in band 1
  // maxNbrParticlesBand0 = maximum number of particles in band 0
  // maxNbrParticlesBand1 = maximum number of particles in band 1
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int kxMomentum, int kyMomentum, char* outputDirectory = 0, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace(const FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace& operator = (const FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // provide the default name of the file for Hilbert space storage 
  //
  // return value = pointer to file name (0 if no default file name exists) 
  virtual char* GetDefaultHilbertSpaceFileName();
  
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
  virtual long GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, unsigned long* singleBandStateDescription, long pos);

  // evaluate all the single band Hilbert spaces
  //
  // maxBandOccupation = maiximum occupation of a single band
  // singleBandTotalKxMax = reference on the array for the maximum total Kx values 
  // singleBandTotalKyMax = reference on the array for the maximum total Ky values
  // singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky)
  // singleBandStates = reference on the array for all the Hilbert space basis states 
  virtual void GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, unsigned long****& singleBandStates);

  // generate all states using the single band hilbert spaces
  //
  virtual void GenerateStatesFromSingleBandHilbertSpaces();
  
};


#endif


