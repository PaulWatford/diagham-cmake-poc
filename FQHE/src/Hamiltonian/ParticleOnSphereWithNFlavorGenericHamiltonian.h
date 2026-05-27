////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of N-flavor Hamiltonian for particles on a sphere            //
//                                                                            //
//                           class author: Sahana Das                         //
//                                                                            //
//                        last modification : 27/05/2026                      //
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



#ifndef PARTICLEONSPHEREWITHNFLAVORGENERICHAMILTONIAN_H
#define PARTICLEONSPHEREWITHNFLAVORGENERICHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithNFlavor.h"
#include "Hamiltonian/AbstractQHEOnSphereWithNFlavorHamiltonian.h"

#include <iostream>

using std::ostream;

class MathematicaOutput;
class AbstractArchitecture;
class QHEParticlePrecalculationOperation;

class ParticleOnSphereWithNFlavorGenericHamiltonian
  : public AbstractQHEOnSphereWithNFlavorHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

protected:

//  double** PseudoPotentials;       // [NbrChannels][LzMax+1]
//  double*** OneBodyPotentials;     // [NbrFlavors][NbrFlavors][LzMax+1]


public:
  ParticleOnSphereWithNFlavorGenericHamiltonian( ParticleOnSphereWithNFlavor* particles, int nbrParticles, int lzmax, double** pseudoPotential, double*** oneBodyPotential, AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  virtual ~ParticleOnSphereWithNFlavorGenericHamiltonian();
  AbstractHamiltonian* Clone ();
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);
  AbstractHilbertSpace* GetHilbertSpace ();

  int GetHilbertSpaceDimension ();

  void ShiftHamiltonian(double shift);

  Complex MatrixElement (RealVector& V1, RealVector& V2);
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  List<Matrix*> LeftInteractionOperators();  
  List<Matrix*> RightInteractionOperators();  

protected:

  virtual void EvaluateInteractionFactors();

};

#endif
