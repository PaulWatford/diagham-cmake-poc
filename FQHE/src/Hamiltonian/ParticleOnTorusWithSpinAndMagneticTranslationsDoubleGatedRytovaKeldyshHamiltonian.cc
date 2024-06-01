////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//   class of hamiltonian associated to spinful particles on a torus with     //
//     double gated Rytova-Keldysh interaction and magnetic translations      //
//                                                                            //
//                        last modification : 06/02/2024                      //
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


#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Architecture/AbstractArchitecture.h"

#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333

// default constructor
//

ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian()
{
  this->ScreeningLength = 0.0;
}


// constructor from pseudopotentials
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// landauLevel = landauLevel to be simulated (GaAs (>=0) or graphene (<0))
// screeningLength = screening length (half the distance between the two screening gates)
// interlayerDistance = interlayer distance a.k.a. the alpha factor in 1 / (q (1 + alpha q))
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

ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int landauLevel, 
																				     double screeningLength, double interlayerDistance, double scalingFactorUpUp, double scalingFactorDownDown, double scalingFactorUpDown,
																				     int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
																				     int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
																				     int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
																				     
																				     double spinFluxUp, double spinFluxDown, 
																				     AbstractArchitecture* architecture, long memory, char* precalculationFileName, double* oneBodyPotentielUpUp, double* oneBodyPotentielDownDown, double* oneBodyPotentielUpDown)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  this->SpinFluxUp = spinFluxUp;
  this->SpinFluxDown = spinFluxDown;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  this->ScreeningLength = screeningLength;
  this->InterlayerDistance = interlayerDistance;
  this->ScalingFactorUpUp = scalingFactorUpUp;
  this->ScalingFactorDownDown = scalingFactorDownDown;
  this->ScalingFactorUpDown = scalingFactorUpDown;
  this->NbrPseudopotentialsUpUp = nbrPseudopotentialsUpUp;
  if (this->NbrPseudopotentialsUpUp > 0)
    {
      this->PseudopotentialsUpUp = new double[this->NbrPseudopotentialsUpUp];
      for (int i = 0; i < this->NbrPseudopotentialsUpUp; ++i)
	this->PseudopotentialsUpUp[i] = pseudopotentialsUpUp[i];
    }
  else
    {
      this->PseudopotentialsUpUp = 0;
    }
  this->NbrPseudopotentialsDownDown = nbrPseudopotentialsDownDown;
  if (this->NbrPseudopotentialsDownDown > 0)
    {
      this->PseudopotentialsDownDown = new double[this->NbrPseudopotentialsDownDown];
      for (int i = 0; i < this->NbrPseudopotentialsDownDown; ++i)
	this->PseudopotentialsDownDown[i] = pseudopotentialsDownDown[i];
    }
  else
    {
      this->PseudopotentialsDownDown = 0;
    }
  this->NbrPseudopotentialsUpDown = nbrPseudopotentialsUpDown;
  if (this->NbrPseudopotentialsUpDown > 0)
    {
      this->PseudopotentialsUpDown = new double[this->NbrPseudopotentialsUpDown];
      for (int i = 0; i < this->NbrPseudopotentialsUpDown; ++i)
	this->PseudopotentialsUpDown[i] = pseudopotentialsUpDown[i];
    }
  else
    {
      this->PseudopotentialsUpDown = 0;
    }
  this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpUp;
  if (this->NbrPseudopotentialsDownDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsDownDown;
  if (this->NbrPseudopotentialsUpDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpDown;
  if (this->MaxNbrPseudopotentials > 0)
    {
      this->LaguerrePolynomials = new Polynomial[this->MaxNbrPseudopotentials];
      for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
        this->LaguerrePolynomials[i] = LaguerrePolynomial(i);
    }
  else
    {
      this->LaguerrePolynomials = 0;
    }
  
  this->LandauLevel = landauLevel;
  if (this->LandauLevel >= 0)
    {
      // simple coulomb interactions
      this->FormFactor = LaguerrePolynomial(this->LandauLevel);
    }
  else
    {
      // coulomb interactions in graphene
      this->FormFactor = 0.5*(LaguerrePolynomial(abs(this->LandauLevel))+LaguerrePolynomial(abs(this->LandauLevel)-1));
    }

  this->OneBodyInteractionFactorsupup = 0;
  if(oneBodyPotentielUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsupup[i] = oneBodyPotentielUpUp[i];
    }
  this->OneBodyInteractionFactorsdowndown = 0;
  if(oneBodyPotentielDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsdowndown[i] = oneBodyPotentielDownDown[i];
    } 
  this->OneBodyInteractionFactorsupdown = 0;
  if(oneBodyPotentielUpDown != 0)
    {
      this->OneBodyInteractionFactorsupdown = new Complex[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsupdown[i] = oneBodyPotentielUpDown[i];
	}
    } 

  this->EvaluateExponentialFactors();
  this->EvaluateInteractionFactors();

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024l)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1l << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1l << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  cout << endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian::~ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian() 
{
}


// get fourier transform of interaction
//
// Q2_half = one half of q² value

double ParticleOnTorusWithSpinAndMagneticTranslationsDoubleGatedRytovaKeldyshHamiltonian::GetVofQ(double Q2_half)
{
  if (Q2_half != 0.0)
    {
      double Q = sqrt(2.0 * Q2_half);
      return (this->FormFactor(Q2_half) * this->FormFactor(Q2_half)) * (tanh (0.5 * this->ScreeningLength * Q)/ Q);
    }
  else
    {
      return (this->FormFactor(Q2_half) * this->FormFactor(Q2_half) * (0.5 * this->ScreeningLength));
    }
}

