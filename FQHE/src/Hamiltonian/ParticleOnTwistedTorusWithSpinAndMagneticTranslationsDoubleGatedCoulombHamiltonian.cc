////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
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


#include "Hamiltonian/ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian.h"
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

ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian::ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian()
{
  this->ScreeningLength = 0.0;
}


// constructor from pseudopotentials
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// ratio = torus aspect ratio (Lx/Ly)
// angle =  angle (in pi units) between the two fundamental cycles of the torus, along (Lx sin theta, Lx cos theta) and (0, Ly)
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

ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian::ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
																				       double ratio, double angle,
																				       double screeningLength, double scalingFactorUpUp, double scalingFactorDownDown, double scalingFactorUpDown,
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
  this->Theta = angle * M_PI;
  this->CosTheta = cos(this->Theta);
  this->SinTheta = sin(this->Theta);
  this->SpinFluxUp = spinFluxUp;
  this->SpinFluxDown = spinFluxDown;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  this->ScreeningLength = screeningLength;
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

ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian::~ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian() 
{
}

// get fourier transform of interaction
//
// Q2_half = one half of q² value

// get fourier transform of the interaction
//
// Q2_half = one half of q² value
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// nonpseudoScaling = rescaling factor for any non-pseudopotential interaction
// return value = Fourrier transform of the interaction

double ParticleOnTwistedTorusWithSpinAndMagneticTranslationsDoubleGatedCoulombHamiltonian::GetVofQ(double q2_half, int nbrPseudopotentials, double* pseudopotentials, double nonpseudoScaling)
{
  double Result = 0.0;
  double Q2 = 2.0 * q2_half;
  if (q2_half != 0.0)
    {
      double Q = sqrt(2.0 * q2_half);
      Result += nonpseudoScaling * (tanh (0.5 * this->ScreeningLength * Q)/ Q);
    }
  else
    {
      Result += nonpseudoScaling * (0.5 * this->ScreeningLength);
    }
  for (int i = 0; i < nbrPseudopotentials; ++i)
    {
      if (pseudopotentials[i] != 0.0)
	{
	  Result += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(Q2);
	}
    }
  return Result * exp(-q2_half);
}

