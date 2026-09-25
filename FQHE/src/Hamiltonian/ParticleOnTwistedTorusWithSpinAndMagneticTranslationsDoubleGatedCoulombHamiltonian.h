////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of hamiltonian associated to spinful particles on a twisted torus   //
//     with double gated coulomb interaction and magnetic translations        //
//                                                                            //
//                        last modification : 15/02/2024                      //
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


#ifndef PARTICLEONTWISTEDTORUSWITHSPINANDMAGNETICTRANSLATIONSDOUBLEGATEDCOULOMBHAMILTONIAN_H
#define PARTICLEONTWISTEDTORUSWITHSPINANDMAGNETICTRANSLATIONSDOUBLEGATEDCOULOMBHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian : public ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonian
{

 protected:

  // screening length (half the distance between the two screening gates)
  double ScreeningLength;

 public:
   
  // default constructor
  // 
  ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian();


  // constructor from pseudopotentials
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the lengths of the two fundamental cycles of the torus L1 / L2
  // angle = angle between the two fundamental cycles of the torus, along (L1 sin, L1 cos) and (0, L2)
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
  ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
									      double ratio, double angle,
									      double screeningLength, double scalingFactorUpUp, double scalingFactorDownDown, double scalingFactorUpDown,
									      int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
									      int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
									      int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
									      double spinFluxUp, double spinFluxDown, 
									      AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
									      double* oneBodyPotentielUpUp = 0, double* oneBodyPotentielDownDown = 0, double* oneBodyPotentielUpDown = 0);
  
  // destructor
  //
  ~ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

 protected:
 
  // get fourier transform of interaction
  //
  // Q2_half = one half of q² value
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // nonpseudoScaling = rescaling factor for any non-pseudopotential interaction
  // return value = Fourrier transform of the interaction
  virtual double GetVofQ(double Q2_half, int nbrPseudopotentials, double* pseudopotentials, double nonpseudoScaling);

};

#endif
