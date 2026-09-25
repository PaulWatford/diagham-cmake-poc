////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//               class of hamiltonian associated to particles                 //
//          on a twisted torus with a generic two body interaction            //
//             and a time reversal symmetric degree of freedom                //
//                                                                            //
//                        last modification : 16/02/2024                      //
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

#ifndef PARTICLEONTWISTEDTORUSWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICGENERICHAMILTONIAN_H
#define PARTICLEONTWISTEDTORUSWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTwistedTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian : public ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonian
{
  
protected:
  
public:
  
  // default constructor    
  //
  ParticleOnTwistedTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian();
  
  // constructor from pseudopotentials
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the lengths of the two fundamental cycles of the torus L1 / L2
  // angle = angle between the two fundamental cycles of the torus, along (L1 sin, L1 cos) and (0, L2)
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
  ParticleOnTwistedTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
											       double ratio, double angle,
											       int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
											       int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
											       int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
											       double spinFluxUp, double spinFluxDown, 
											       AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
											       double* oneBodyPotentielUpUp = 0, double* oneBodyPotentielDownDown = 0, double* oneBodyPotentielUpDown = 0);
  
  // destructor
  //
  ~ParticleOnTwistedTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian();
  
  protected:
  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

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
  virtual Complex ComplexEvaluateInteractionCoefficientUpUp(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
							    double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);
  
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
  virtual Complex ComplexEvaluateInteractionCoefficientDownDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
								double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);
  
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
  virtual Complex ComplexEvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
							      double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4); 

  // get all the indices that should appear in the annihilation/creation operators
  //
  virtual void GetIndices();
  
};

#endif
