////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of a two body interaction projected onto three bands with       //
//      spin-like degree of freedom (requiring only U(1) conservation)        //
//        from an ASCII file providing the two body matrix elements           //
//                                                                            //
//                        last modification : 16/11/2023                      //
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


#ifndef PARTICLEONLATTICEFROMFILEINTERACTIONTHREEBANDWITHSPINHAMILTONIAN_H
#define PARTICLEONLATTICEFROMFILEINTERACTIONTHREEBANDWITHSPINHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian : public ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian
{

 protected:
 
 public:

  // default constructor
  //
  ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // matrixElementsInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // interactionRescalingFactor = global rescaling factor for the two-body interaction term
  // additionalSpinFlag = include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
								   char* matrixElementsInteractionFile,
								   Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag,
								   double interactionRescalingFactor, 
								   bool additionalSpinFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian();
  

 protected:
 
  // convert spin and band indices into a linearized index
  //
  // spinValue = spin value (+1 or -1)
  // bandIndex = band index
  // return value = linearized index
  virtual int GetLinearizedSpinBandIndex(int spinValue, int bandIndex);

  // test if the sum of spin projection (s_1+s2==s_3+s_4) is conserved from the linearized indices of the internal degrees of freedom
  //
  // sigma1 = first linearized index
  // sigma2 = second linearized index
  // sigma3 = third linearized index
  // sigma4 = fourth linearized index
  // return value = true if the spin projection is conserved
  virtual bool TestSpinConservation(int sigma1, int sigma2, int sigma3, int sigma4);

  // test if the sum of spin projection (s_1+s2==s_3+s_4) is conserved from the linearized indices of the internal degrees of freedom (assuming both spin and valley)
  //
  // sigma1 = first linearized index
  // sigma2 = second linearized index
  // sigma3 = third linearized index
  // sigma4 = fourth linearized index
  // return value = true if the spin projection is conserved
  virtual bool TestSpinfulValleyConservation(int sigma1, int sigma2, int sigma3, int sigma4);

};

// convert spin and band indices into a linearized index
//
// spinValue = spin value (+1 or -1)
// bandIndex = band index
// return value = linearized index

inline int ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian::GetLinearizedSpinBandIndex(int spinValue, int bandIndex)
{
  return ((3 * ((spinValue + 1) >> 1))+ bandIndex);
}

// test if the sum of spin projection (s_1+s2==s_3+s_4) is conserved from the linearized indices of the internal degrees of freedom
//
// sigma1 = first linearized index
// sigma2 = second linearized index
// sigma3 = third linearized index
// sigma4 = fourth linearized index
// return value = true if the spin projection is conserved

inline bool ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian::TestSpinConservation(int sigma1, int sigma2, int sigma3, int sigma4)
{
  return ((((sigma1 / 3) & 1) + ((sigma2 / 3) & 1)) == (((sigma3 / 3) & 1) + ((sigma4 / 3) & 1)));
}


// test if the sum of spin projection (s_1+s2==s_3+s_4) is conserved from the linearized indices of the internal degrees of freedom (assuming both spin and valley)
//
// sigma1 = first linearized index
// sigma2 = second linearized index
// sigma3 = third linearized index
// sigma4 = fourth linearized index
// return value = true if the spin projection is conserved

inline bool ParticleOnLatticeFromFileInteractionThreeBandWithSpinHamiltonian::TestSpinfulValleyConservation(int sigma1, int sigma2, int sigma3, int sigma4)
{
  return  (((sigma1 / 6) == (sigma3 / 6)) && ((sigma2 / 6) == (sigma4 / 6)));
}


#endif
