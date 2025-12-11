////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class density operator (including non-diagonal terms)       //
//              for particle with any type of internal degree of freedom      //
//                                                                            //
//                        last modification : 13/11/2023                      //
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


#include "config.h"
#include "Operator/ParticleOnSquareLatticeWithGenericSpinBandDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default data
//
// particle = hilbert space associated to the particles
// m1 = momentum index of the creation operator
// sigma1 = internal degree of freedom label of the creation operator
// m2 = momentum index of the annihilation operator
// sigma2 = internal degree of freedom label of the annihilation operator

ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::ParticleOnSquareLatticeWithGenericSpinBandDensityOperator(ParticleOnSphereWithSpin* particle, int m1, int sigma1, int m2, int sigma2)
{
  this->Particle = particle;
  this->Momentum1 = m1;
  this->Sigma1 = sigma1;
  this->Momentum2 = m2;
  this->Sigma2 = sigma2;
}

// destructor
//

ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::~ParticleOnSquareLatticeWithGenericSpinBandDensityOperator()
{
}
  
// clone operator without duplicating data
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      int Index = this->Particle->AdsigmaAsigma(i, this->Momentum1, this->Sigma1, this->Momentum2, this->Sigma2, Coefficient);
      if (Index != FullDim)
	{
	  Element += V1[Index] * V2[i] * Coefficient;
	}
    }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
											   int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      int Index = this->Particle->AdsigmaAsigma(i, this->Momentum1, this->Sigma1, this->Momentum2, this->Sigma2, Coefficient);
      if (Index != Dim)
	{
	  vDestination[Index] += vSource[i] * Coefficient;
	}
    }
  return vDestination;
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      int Index = this->Particle->AdsigmaAsigma(i, this->Momentum1, this->Sigma1, this->Momentum2, this->Sigma2, Coefficient);
      if (Index != FullDim)
	{
	  Element += Conj(V1[Index]) * V2[i] * Coefficient;
	}
    }
  return Element;
}
  
// evaluate part of the matrix element without complex conjugate for the left vector, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::ConjugatePartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      int Index = this->Particle->AdsigmaAsigma(i, this->Momentum1, this->Sigma1, this->Momentum2, this->Sigma2, Coefficient);
      if (Index != FullDim)
	{
	  Element += V1[Index] * V2[i] * Coefficient;
	}
    }
  return Element;
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSquareLatticeWithGenericSpinBandDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      int Index = this->Particle->AdsigmaAsigma(i, this->Momentum1, this->Sigma1, this->Momentum2, this->Sigma2, Coefficient);
      if (Index != Dim)
	{
	  vDestination[Index] += vSource[i] * Coefficient;      
	}
    }
  return vDestination;
}

