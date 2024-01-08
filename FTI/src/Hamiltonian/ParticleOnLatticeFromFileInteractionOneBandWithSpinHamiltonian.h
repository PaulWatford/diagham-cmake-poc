////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of a two body interaction projected onto a single band with     //
//      spin-like degree of freedom (requiring only U(1) conservation)        //
//        from an ASCII file providing the two body matrix elements           //
//                                                                            //
//                        last modification : 30/12/2023                      //
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


#ifndef PARTICLEONLATTICEFROMFILEINTERACTIONONEBANDWITHSPINHAMILTONIAN_H
#define PARTICLEONLATTICEFROMFILEINTERACTIONONEBANDWITHSPINHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian : public ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian
{

 protected:
 
 public:

  // default constructor
  //
  ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian();

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
  ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
								 char* matrixElementsInteractionFile,
								 Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag,
								 double interactionRescalingFactor, 
								 bool additionalSpinFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian();
  

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

  // process the matrix elements from the ascii file
  //
  // arraySigma1 = reference on the array containing the indices of the internal degree freedom for the operator 1
  // arraySigma2 = reference on the array containing the indices of the internal degree freedom for the operator 2
  // arraySigma3 = reference on the array containing the indices of the internal degree freedom for the operator 3
  // arraySigma4 = reference on the array containing the indices of the internal degree freedom for the operator 4
  // arraySpinIndex1 = reference on the array containing the indices of the spin degree freedom for the operator 1
  // arraySpinIndex2 = reference on the array containing the indices of the spin degree freedom for the operator 2
  // arraySpinIndex3 = reference on the array containing the indices of the spin degree freedom for the operator 3
  // arraySpinIndex4 = reference on the array containing the indices of the spin degree freedom for the operator 4
  // arrayKx1 = reference on the array containing the momentum along x for the operator 1
  // arrayKy1 = reference on the array containing the momentum along y for the operator 1
  // arrayKx2 = reference on the array containing the momentum along x for the operator 2
  // arrayKy2 = reference on the array containing the momentum along y for the operator 2
  // arrayKx3 = reference on the array containing the momentum along x for the operator 3
  // arrayKy3 = reference on the array containing the momentum along y for the operator 3
  // arrayKx4 = reference on the array containing the momentum along x for the operator 4
  // arrayKy4 = reference on the array containing the momentum along y for the operator 4
  // arrayMatrixElements = reference on the array containing the matrix elements
  // return value = number of entries
  virtual int ProcessTwoBodyMatrixElements(int*& arraySigma1, int*& arraySigma2, int*& arraySigma3, int*& arraySigma4, int*& arraySpinIndex1, int*& arraySpinIndex2, int*& arraySpinIndex3, int*& arraySpinIndex4, int*& arrayKx1, int*& arrayKy1, int*& arrayKx2, int*& arrayKy2, int*& arrayKx3, int*& arrayKy3, int*& arrayKx4, int*& arrayKy4, Complex*& arrayMatrixElements);

};

// convert spin and band indices into a linearized index
//
// spinValue = spin value (+1 or -1)
// bandIndex = band index
// return value = linearized index

inline int ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian::GetLinearizedSpinBandIndex(int spinValue, int bandIndex)
{
  return ((spinValue + 1) >> 1);
}

// test if the sum of spin projection (s_1+s2==s_3+s_4) is conserved from the linearized indices of the internal degrees of freedom
//
// sigma1 = first linearized index
// sigma2 = second linearized index
// sigma3 = third linearized index
// sigma4 = fourth linearized index
// return value = true if the spin projection is conserved

inline bool ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian::TestSpinConservation(int sigma1, int sigma2, int sigma3, int sigma4)
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

inline bool ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian::TestSpinfulValleyConservation(int sigma1, int sigma2, int sigma3, int sigma4)
{
  return  (((sigma1 >> 1) == (sigma3 >> 1)) && ((sigma2 >> 1) == (sigma4 >> 1)));
}


#endif
