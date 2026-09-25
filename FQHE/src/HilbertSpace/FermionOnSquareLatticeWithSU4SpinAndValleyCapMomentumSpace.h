////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//           class of fermions on a square lattice with SU(4) spin            //
//   in momentum space with a cap on the number of particles in each valley   //
//                                                                            //
//                        last modification : 17/11/2023                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU4SPINANDVALLEYCAPMOMENTUMSPACE_H
#define FERMIONONSQUARELATTICEWITHSU4SPINANDVALLEYCAPMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"


#include <iostream>



class FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace : public FermionOnSquareLatticeWithSU4SpinMomentumSpace
{

 protected:

  // minimum number of particles in valley plus (or band 0)
  int MinNbrParticlesPlus;
  // minimum number of particles in valley minus (or band 1)
  int MinNbrParticlesMinus;
  
  // maximum number of particles in valley plus (or band 0)
  int MaxNbrParticlesPlus;
  // maximum number of particles in valley minus (or band 1)
  int MaxNbrParticlesMinus;

 public:

  // default constructor
  //
  FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // maxNbrParticlesPlus = maximum number of particles in valley plus
  // maxNbrParticlesMinus = maximum number of particles in valley minus
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesPlus, int maxNbrParticlesMinus, int kxMomentum, int kyMomentum, char* outputDirectory = 0, unsigned long memory = 10000000);
  
  // constructor when preserving only spin
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // maxNbrParticlesPlus = maximum number of particles in valley plus
  // maxNbrParticlesMinus = maximum number of particles in valley minus
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // totalSpin = twice the total spin value
  // outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesPlus, int maxNbrParticlesMinus, int kxMomentum, int kyMomentum, int totalSpin, char* outputDirectory = 0, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace(const FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace& operator = (const FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // provide the default name of the file for Hilbert space storage 
  //
  // return value = pointer to file name (0 if no default file name exists) 
  virtual char* GetDefaultHilbertSpaceFileName();
  
 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // maxNbrParticlesPlus = current maximum number of particles in valley plus
  // maxNbrParticlesMinus = current maximum number of particles in valley minus
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesPlus, int maxNbrParticlesMinus);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrFermionsUp = current number of fermions with a spin up
  // maxNbrParticlesPlus = current maximum number of particles in valley plus
  // maxNbrParticlesMinus = current maximum number of particles in valley minus
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrFermionsUp, int maxNbrParticlesPlus, int maxNbrParticlesMinus);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // maxNbrParticlesPlus = current maximum number of particles in valley plus
  // maxNbrParticlesMinus = current maximum number of particles in valley minus
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int maxNbrParticlesPlus, int maxNbrParticlesMinus);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // nbrFermionsUp = current number of fermions with a spin up
  // maxNbrParticlesPlus = current maximum number of particles in valley plus
  // maxNbrParticlesMinus = current maximum number of particles in valley minus
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int nbrFermionsUp, int maxNbrParticlesPlus, int maxNbrParticlesMinus);

  // generate all states using the single band hilbert spaces
  //
  virtual void GenerateStatesFromSingleBandHilbertSpaces();

};


#endif


