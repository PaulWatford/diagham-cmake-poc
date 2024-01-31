////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//   class of hamiltonian associated to spinful particles on a torus with     //
//         double gated coulomb interaction, time reversal symmetry           //
//                           and magnetic-like translations                   //
//                                                                            //
//                        last modification : 30/01/2024                      //
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


#ifndef PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICDOUBLEGATEDCOULOMBHAMILTONIAN_H
#define PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICDOUBLEGATEDCOULOMBHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricDoubleGatedCoulombHamiltonian :public ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian
{

 protected:

  // screening length (half the distance between the two screening gates)
  double ScreeningLength;

 public:
      
  // default constructor
  //
  ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricDoubleGatedCoulombHamiltonian();
      
  // constructor from pseudopotentials
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // screeningLength = screening length (half the distance between the two screening gates)
  // scalingFactorUpUp = global rescaling factor for the up-up interaction
  // scalingFactorDownDown = global rescaling factor for the down-down interaction
  // scalingFactorUpDown = global rescaling factor for the up-down interaction
  // nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
  // pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
  // nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
  // pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
  // nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
  // pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
  // spinFluxUp = additional inserted flux for spin up
  // spinFluxDown = additional inserted flux for spin down
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricDoubleGatedCoulombHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
												   double ratio,
												   double screeningLength, double scalingFactorUpUp, double scalingFactorDownDown, double scalingFactorUpDown,
												   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
												   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
												   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
												   double spinFluxUp, double spinFluxDown, 
												   AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
												   double* oneBodyPotentielUpUp = 0, double* oneBodyPotentielDownDown = 0, double* oneBodyPotentielUpDown = 0);
  
  // destructor
  //
  ~ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricDoubleGatedCoulombHamiltonian();

  
  protected:
 
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // spinFluxM1 = additional inserted flux for m1
  // spinFluxM2 = additional inserted flux for m2
  // spinFluxM3 = additional inserted flux for m3
  // spinFluxM4 = additional inserted flux for m4
  // nonpseudoScaling = rescaling factor for any non-pseudopotential interaction
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
						double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4, double nonpseudoScaling);
  
  
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // spinFluxM1 = additional inserted flux for m1
  // spinFluxM2 = additional inserted flux for m2
  // spinFluxM3 = additional inserted flux for m3
  // spinFluxM4 = additional inserted flux for m4
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
						      double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4); 
  
  // get fourier transform of interaction
  //
  // Q2_half = one half of q² value  
  virtual double GetVofQ(double Q2_half);

  };

#endif
