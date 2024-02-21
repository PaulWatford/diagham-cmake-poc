////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
// class density-density operator (normal ordered and including non-diagonal  //
//      terms) for particle with any type of internal degree of freedom       //
//                                                                            //
//                        last modification : 18/02/2024                      //
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
#include "Operator/ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default data
//
// particle = hilbert space associated to the particles
// m1 = momentum index of the leftmost creation operator
// sigma1 = internal degree of freedom label of the leftmost creation operator
// m2 = momentum index of the rightmost creation operator
// sigma2 = internal degree of freedom label of the rightmost creation operator
// m3 = momentum index of the leftmost annihilation operator
// sigma3 = internal degree of freedom label of the leftmost annihilation operator
// m4 = momentum index of the rightmost annihilation operator
// sigma4 = internal degree of freedom label of the rightmost annihilation operator

ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator(ParticleOnSphereWithSpin* particle, int m1, int sigma1, int m2, int sigma2, int m3, int sigma3, int m4, int sigma4)
{
  this->Particle = particle;
  this->Momentum1 = m1;
  this->Sigma1 = sigma1;
  this->Momentum2 = m2;
  this->Sigma2 = sigma2;
  this->Momentum3 = m3;
  this->Sigma3 = sigma3;
  this->Momentum4 = m4;
  this->Sigma4 = sigma4;
}

// destructor
//

ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::~ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator()
{
}
  
// clone operator without duplicating data
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      Coefficient = this->Particle->AsigmaAsigma(i, this->Momentum3, this->Momentum4, this->Sigma3, this->Sigma4);
      if (Coefficient != 0.0)
	{
	  int Index = this->Particle->AdsigmaAdsigma(this->Momentum1, this->Momentum2, this->Sigma1, this->Sigma2, Coefficient);
	  if (Index != FullDim)
	    {
	      Element += V1[Index] * V2[i] * Coefficient;
	    }
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

RealVector& ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
												  int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Coefficient = this->Particle->AsigmaAsigma(i, this->Momentum3, this->Momentum4, this->Sigma3, this->Sigma4);
      if (Coefficient != 0.0)
	{
	  int Index = this->Particle->AdsigmaAdsigma(this->Momentum1, this->Momentum2, this->Sigma1, this->Sigma2, Coefficient);
	  if (Index != FullDim)
	    {
	      vDestination[Index] += vSource[i] * Coefficient;
	    }
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

Complex ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      Coefficient = this->Particle->AsigmaAsigma(i, this->Momentum3, this->Momentum4, this->Sigma3, this->Sigma4);
      if (Coefficient != 0.0)
	{
	  int Index = this->Particle->AdsigmaAdsigma(this->Momentum1, this->Momentum2, this->Sigma1, this->Sigma2, Coefficient);
	  if (Index != FullDim)
	    {
	      Element += Conj(V1[Index]) * V2[i] * Coefficient;
	    }
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

ComplexVector& ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Coefficient = this->Particle->AsigmaAsigma(i, this->Momentum3, this->Momentum4, this->Sigma3, this->Sigma4);
      if (Coefficient != 0.0)
	{
	  int Index = this->Particle->AdsigmaAdsigma(this->Momentum1, this->Momentum2, this->Sigma1, this->Sigma2, Coefficient);
	  if (Index != FullDim)
	    {
	      vDestination[Index] += vSource[i] * Coefficient;      
	    }
	}
    }
  return vDestination;
}

