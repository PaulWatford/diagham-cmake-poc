////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//        class of a two body interaction projected onto a single band        //
//        from an ASCII file providing the three body matrix elements         //
//                                                                            //
//                        last modification : 07/05/2024                      //
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


#ifndef PARTICLEONLATTICEFROMFILEINTERACTIONONEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEFROMFILEINTERACTIONONEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullTwoBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian : public ParticleOnLatticeQuantumSpinHallFullTwoBandNBodyHamiltonian
{

 protected:
  
  // pointer to the tight binding model
  Abstract2DTightBindingModel* TightBindingModel;

  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;

  // use flat band model
  bool FlatBand;
  
  // number of internal indices (i.e number of entries for indices 0, 1, 2 and 3 of InteractionFactorsSigma)
  int NbrInternalIndices;

  // global rescaling factor for the three-body interaction term
  double ThreeBodyInteractionRescalingFactor;

  // name of the ASCII file containing the matrix element for the generic three body interaction term
  char* MatrixElementsThreeBodyInteractionFile;

  // global rescaling factor for the two-body interaction term
  double TwoBodyInteractionRescalingFactor;

  // name of the ASCII file containing the matrix element for the generic two body interaction term
  char* MatrixElementsTwoBodyInteractionFile;


 public:

  // default constructor
  //
  ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // matrixElementsThreeBodyInteractionFile = name of the ASCII file containing the matrix element for the generic three body interaction term
  // threeBodyInteractionRescalingFactor = global rescaling factor for the three-body interaction term
  // matrixElementsTwoBodyInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
  // twoBodyInteractionRescalingFactor = global rescaling factor for the two-body interaction term
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // spinFlag = include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,	
								  char* matrixElementsThreeBodyInteractionFile, double threeBodyInteractionRescalingFactor,
								  char* matrixElementsTwoBodyInteractionFile,  double twoBodyInteractionRescalingFactor, 
								  Abstract2DTightBindingModel* tightBindingModel ,bool flatBandFlag,
								  bool spinFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate the three-body interaction factors 
  //   
  virtual void EvaluateThreeBodyInteractionFactors();
  
  // evaluate  the two-body interaction factors
  //   
  virtual void EvaluateTwoBodyInteractionFactors();
  
  // evaluate all one-body factors
  //   
  virtual void EvaluateOneBodyInteractionFactors();

  // process the three-body matrix elements from the ascii file
  //
  // arraySigma1 = reference on the array containing the indices of the internal degree freedom for the operator 1
  // arraySigma2 = reference on the array containing the indices of the internal degree freedom for the operator 2
  // arraySigma3 = reference on the array containing the indices of the internal degree freedom for the operator 3
  // arraySigma4 = reference on the array containing the indices of the internal degree freedom for the operator 4
  // arraySigma5 = reference on the array containing the indices of the internal degree freedom for the operator 5
  // arraySigma6 = reference on the array containing the indices of the internal degree freedom for the operator 6
  // arrayKx1 = reference on the array containing the momentum along x for the operator 1
  // arrayKy1 = reference on the array containing the momentum along y for the operator 1
  // arrayKx2 = reference on the array containing the momentum along x for the operator 2
  // arrayKy2 = reference on the array containing the momentum along y for the operator 2
  // arrayKx3 = reference on the array containing the momentum along x for the operator 3
  // arrayKy3 = reference on the array containing the momentum along y for the operator 3
  // arrayKx4 = reference on the array containing the momentum along x for the operator 4
  // arrayKy4 = reference on the array containing the momentum along y for the operator 4
  // arrayKx5 = reference on the array containing the momentum along x for the operator 5
  // arrayKy5 = reference on the array containing the momentum along y for the operator 5
  // arrayKx6 = reference on the array containing the momentum along x for the operator 6
  // arrayKy6 = reference on the array containing the momentum along y for the operator 6
  // arrayMatrixElements = reference on the array containing the matrix elements
  // return value = number of entries
  virtual int ProcessThreeBodyMatrixElements(int*& arraySigma1, int*& arraySigma2, int*& arraySigma3, int*& arraySigma4, int*& arraySigma5, int*& arraySigma6,
					     int*& arrayKx1, int*& arrayKy1, int*& arrayKx2, int*& arrayKy2, int*& arrayKx3, int*& arrayKy3,
					     int*& arrayKx4, int*& arrayKy4, int*& arrayKx5, int*& arrayKy5, int*& arrayKx6, int*& arrayKy6, Complex*& arrayMatrixElements);

  // process the matrix elements from the ascii file
  //
  // arraySigma1 = reference on the array containing the indices of the internal degree freedom for the operator 1
  // arraySigma2 = reference on the array containing the indices of the internal degree freedom for the operator 2
  // arraySigma3 = reference on the array containing the indices of the internal degree freedom for the operator 3
  // arraySigma4 = reference on the array containing the indices of the internal degree freedom for the operator 4
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
  virtual int ProcessTwoBodyMatrixElements(int*& arraySigma1, int*& arraySigma2, int*& arraySigma3, int*& arraySigma4, int*& arrayKx1, int*& arrayKy1, int*& arrayKx2, int*& arrayKy2, int*& arrayKx3, int*& arrayKy3, int*& arrayKx4, int*& arrayKy4, Complex*& arrayMatrixElements);

  // test if the sum of spin projection is conserved from the linearized indices of the internal degrees of freedom (assuming both spin and valley)
  //
  // sigma1 = first linearized index
  // sigma2 = second linearized index
  // sigma3 = third linearized index
  // sigma4 = fourth linearized index
  // return value = true if the spin projection is conserved
  virtual bool TestSpinfulValleyConservation(int sigma1, int sigma2, int sigma3, int sigma4);
  
  // initialize the n-body interaction terms
  //
  virtual void InitializeNBodyInteraction();

};


// test if the sum of spin projection is conserved from the linearized indices of the internal degrees of freedom (assuming both spin and valley)
//
// sigma1 = first linearized index
// sigma2 = second linearized index
// sigma3 = third linearized index
// sigma4 = fourth linearized index
// return value = true if the spin projection is conserved

inline bool ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::TestSpinfulValleyConservation(int sigma1, int sigma2, int sigma3, int sigma4)
{
  return  (((sigma1 / 3) == (sigma3 / 3)) && ((sigma2 / 3) == (sigma4 / 3)));
}


#endif
