////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of pair-hopping hamiltonian at momentum points 0 or pi          //
//                                                                            //
//                        last modification : 01/06/2022                      //
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


#include "Hamiltonian/ParticleOnTorusPairHoppingRealHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Polynomial/SpecialPolynomial.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnTorusPairHoppingRealHamiltonian::ParticleOnTorusPairHoppingRealHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusPairHoppingRealHamiltonian::ParticleOnTorusPairHoppingRealHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
										     int nbrParticles, int nbrSites, int xMomentum,
										     AbstractArchitecture* architecture, long memory, 
										     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = nbrSites - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = nbrSites;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyTermFlag = false;
  this->OneBodyInteractionFactors = 0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  this->HamiltonianShift = 0.0;
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
      else
	{
	  if (this->Architecture->HasAutoLoadBalancing())
	    {
	      this->FastMultiplicationMemory(0l);
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusPairHoppingRealHamiltonian::~ParticleOnTorusPairHoppingRealHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusPairHoppingRealHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusPairHoppingRealHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusPairHoppingRealHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();
  this->InteractionFactors = new double* [this->NbrSectorSums];
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      double TmpPHEnergyShift = 0.0;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	      TmpPHEnergyShift += (this->EvaluateInteractionCoefficient(m1, m2, m1, m2)
				   + this->EvaluateInteractionCoefficient(m2, m1, m2, m1)
				   - this->EvaluateInteractionCoefficient(m1, m2, m2, m1)
				   - this->EvaluateInteractionCoefficient(m2, m1, m1, m2));
	    }
	}
      cout << "ph energy shift = " << TmpPHEnergyShift << endl;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					   - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					   - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors[i][Index] = TmpCoefficient;
		      TotalNbrNonZeroInteractionFactors++;
		    }
		  else
		    {
		      this->InteractionFactors[i][Index] = 0.0;
		    }
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					   + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (m3 == m4)
		    TmpCoefficient *= 0.5;
		  if (m1 == m2)
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					   + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (m3 == m4)
		    TmpCoefficient *= 0.5;
		  if (m1 == m2)
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors[i][Index] = TmpCoefficient;
		      TotalNbrNonZeroInteractionFactors++;
		    }
		  else
		    {
		      this->InteractionFactors[i][Index] = 0.0;
		    }
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "nbr non-zero interaction = " << TotalNbrNonZeroInteractionFactors << endl;
  cout << "====================================" << endl;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusPairHoppingRealHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  if (m1 > m2)
    {
      m2 += this->MaxMomentum;
    }
  if (m3 > m4)
    {
      m4 += this->MaxMomentum;
    }
  // if ((m1 != m3) && (m1 != m4) && (m2 != m3) && (m2 != m4))
  //   cout << "m1=" << m1 << " " << "m2=" << m2 << " " << "m3=" << m3 << " " << "m4=" << m4 << " " << endl;
  
  if ((((m3 + 1) == m4) && ((m1 + 3) == m2) && (((m1 + 1) == m3) || ((m1 + 1) == (m3 + this->MaxMomentum))))
      || (((m3 + 3) == m4) && ((m1 + 1) == m2) && (((m3 + 1) == m1) || ((m3 + 1) == (m1 + this->MaxMomentum)))))
    {
      //      cout << "check" << endl;
      return 1.0;
    }
  else
    {
      return 0.0;
    }
}

