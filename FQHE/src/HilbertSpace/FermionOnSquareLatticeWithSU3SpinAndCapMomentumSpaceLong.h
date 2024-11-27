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
//                 with a cap on the number of particles per band             //
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


#ifndef FERMIONONSQUARELATTICEWITHSU3SPINANDCAPMOMENTUMSPACELONG_H
#define FERMIONONSQUARELATTICEWITHSU3SPINANDCAPMOMENTUMSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong.h"

#include <iostream>



class FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong : public FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong
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
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong();
  
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
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

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
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong(const FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong& operator = (const FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
									 int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture = 0);
  
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
  virtual long GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, ULONGLONG* singleBandStateDescription, long pos);

  // evaluate all the single band Hilbert spaces
  //
  // maxBandOccupation = maiximum occupation of a single band
  // singleBandTotalKxMax = reference on the array for the maximum total Kx values 
  // singleBandTotalKyMax = reference on the array for the maximum total Ky values
  // singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky)
  // singleBandStates = reference on the array for all the Hilbert space basis states 
  virtual void GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, ULONGLONG****& singleBandStates);

  // generate all states using the single band hilbert spaces
  //
  virtual void GenerateStatesFromSingleBandHilbertSpaces();

};


#endif


