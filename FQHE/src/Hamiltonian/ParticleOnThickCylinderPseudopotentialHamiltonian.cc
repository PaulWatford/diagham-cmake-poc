////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of hamiltonian associated to particles on a           //
//                   thick cylinder with pseudopotential interaction          //
//                                                                            //
//                      class author: Andreas Feuerpfeil                      //
//                                                                            //
//                        last modification : 26/06/2026                      //
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


#include "Hamiltonian/ParticleOnThickCylinderPseudopotentialHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Polynomial/SpecialPolynomial.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = array containing pseudopotential values
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnThickCylinderPseudopotentialHamiltonian::ParticleOnThickCylinderPseudopotentialHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
												     int nbrPseudopotentials, double* pseudopotentials, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = 0.0;
  this->InvRatio = 0.0;

  this->NbrPseudopotentials = nbrPseudopotentials;
  this->Pseudopotentials = new double[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->Pseudopotentials[i] = pseudopotentials[i];
  this->LaguerrePolynomials =new Polynomial[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->Architecture = architecture;
  this->Confinement = 0.0;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;


  this->OneBodyInteractionFactors = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
   { 
       this->OneBodyInteractionFactors[i] = 0.0;
    }

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1 << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1 << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
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

ParticleOnThickCylinderPseudopotentialHamiltonian::~ParticleOnThickCylinderPseudopotentialHamiltonian() 
{
}


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnThickCylinderPseudopotentialHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Length = (double) this->NbrLzValue;
  double kappa = 2.0 / Length;
  double Omega = 2.0 * M_PI / kappa;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;
  double error;
  
  double Coefficient = 0.0;

  Coefficient = this->ThickPseudopotentialMatrixElement(Xm1-Xm4, Xm1-Xm3, this->NbrPseudopotentials, this->Pseudopotentials, this->LaguerrePolynomials, error);

  if (fabs(error) > 1e-6)
    {
      cout << "Warning: large error in matrix elements! " ;
    }

  return (Coefficient/(2.0 * Omega));
}

// evaluate the matrix element of V(q)= sum_i (-1)^i * pseudopotentials[i] * q^(2i) for a thick cylinder
//

double ParticleOnThickCylinderPseudopotentialHamiltonian::ThickPseudopotentialMatrixElement(
											    double q1, double q2, int nbrPseudopotentials,
											    const double* pseudopotentials, Polynomial* laguerrePolynomials, double& error)
{
  (void) q2; 
  (void) laguerrePolynomials;
  error = 0.0;

  double neg_q1_sq = -(q1 * q1);

  // Horner's method: sum_i (-1)^i * pseudopotentials[i] * q1^(2i)
  // Equivalent to evaluating polynomial in x = -q1^2 with coefficients pseudopotentials[i]
  double sum = pseudopotentials[nbrPseudopotentials - 1];
  for (int i = nbrPseudopotentials - 2; i >= 0; --i) {
    sum = sum * neg_q1_sq + pseudopotentials[i];
  }
  return sum;
}
