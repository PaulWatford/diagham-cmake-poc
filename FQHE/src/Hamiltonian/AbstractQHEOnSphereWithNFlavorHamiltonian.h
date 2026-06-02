////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of generic N-flavor Hamiltonian for particles on a           //
//                                    sphere                                  //
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


#ifndef ABSTRACTQHEONSPHEREWITHNFLAVORHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREWITHNFLAVORHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithNFlavor.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>

using std::ostream;



class AbstractQHEOnSphereWithNFlavorHamiltonian: public AbstractQHEOnSphereHamiltonian
{

protected:

  // ============================================================
  // Flavor structure
  // ============================================================

  int NbrFlavors;      // N
  int NbrChannels;     // N(N+1)/2

  inline int Channel(int s1, int s2) const
  {
    /*if (this->NbrFlavors == 2)
    {
      if (s1 == 0 && s2 == 0) return 0;   // up-up
      if (s1 == 1 && s2 == 1) return 1;   // down-down
      return 2;                           // up-down
    }*/
    if (s1 > s2)
      std::swap(s1, s2);
    return s1 * NbrFlavors - (s1*(s1-1))/2 + (s2 - s1);
  }

  // ============================================================
  // Orbital sector grouping
  // ============================================================

 
  int NbrIntraSectorSums;
  int* NbrIntraSectorIndicesPerSum;
  int** IntraSectorIndicesPerSum;

  int NbrInterSectorSums;
  int* NbrInterSectorIndicesPerSum;
  int** InterSectorIndicesPerSum;

  double*** InteractionFactorsIntra;
  double*** InteractionFactorsInter;

  double*** OneBodyInteractionFactors;
  double*** OneBodyInteractionFactorsIntra;
  double*** OneBodyInteractionFactorsInter;

  // ============================================================
  // Alternative fast-multiplication indexing
  // ============================================================

  int* M1IntraValue;
  int* M2IntraValue;
  int NbrM12IntraIndices;
  int** M3IntraValues;
  int* NbrM3IntraValues;

  int* M1InterValue;
  int* M2InterValue;
  int NbrM12InterIndices;
  int** M3InterValues;
  int* NbrM3InterValues;

  // ============================================================
  // Pseudopotentials
  double** PseudoPotentials;       // [NbrChannels][LzMax+1]
  double*** OneBodyPotentials;     // [NbrFlavors][NbrFlavors][LzMax+1]


public:

  AbstractQHEOnSphereWithNFlavorHamiltonian();

  // virtual destructor
  //
  virtual ~AbstractQHEOnSphereWithNFlavorHamiltonian();

  virtual bool IsHermitian();
  virtual bool IsConjugate();

  virtual RealVector& LowLevelAddMultiply( RealVector& vSource, RealVector& vDestination, 
					   int firstComponent, int nbrComponent);


  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations,
						  int nbrVectors, int firstComponent, int nbrComponent);


protected:

  virtual void EvaluateInteractionFactors() = 0;
  virtual void EnableFastMultiplication();

  RealVector* LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations,
							     int nbrVectors, int firstComponent, int nbrComponent);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  virtual void PartialEnableFastMultiplication( int firstComponent, int lastComponent);

  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithNFlavor* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination);
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithNFlavor* particles, int index, int* indexArray, double* coefficientArray, long& position);
};



inline void AbstractQHEOnSphereWithNFlavorHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(
  ParticleOnSphereWithNFlavor* particles,
  int firstComponent, int lastComponent, int step,
  RealVector& vSource,
  RealVector& vDestination)
{

int Dim = particles->GetHilbertSpaceDimension();
double Coefficient;
int Index;

for (int i = firstComponent; i < lastComponent; i += step)
{
  double Source = vSource[i];

  /* Diagonal terms */
  if (this->OneBodyPotentials != 0)
  {
    double TmpDiagonal = 0.0;

    for (int m = 0; m <= this->LzMax; ++m)
      for (int a = 0; a < this->NbrFlavors; ++a)
        TmpDiagonal += this->OneBodyPotentials[a][a][m]
                       * particles->AdsigmaAsigma(i, m, a);

    vDestination[i] += (this->HamiltonianShift + TmpDiagonal) * Source;
  }
  else
    vDestination[i] += this->HamiltonianShift * Source;

  /* Off-diagonal tunneling */
  if (this->OneBodyPotentials != 0)
    for (int m = 0; m <= this->LzMax; ++m)
      for (int a = 0; a < this->NbrFlavors; ++a)
        for (int b = 0; b < this->NbrFlavors; ++b)
          if (a != b)
          {
            Index = particles->AdsigmaAsigma(i, m, m, a, b, Coefficient);
            if (Index < Dim)
              vDestination[Index] += Coefficient
                  * this->OneBodyPotentials[a][b][m]
                  * Source;
          }
}
}



inline void AbstractQHEOnSphereWithNFlavorHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(
  ParticleOnSphereWithNFlavor* particles,
  int index,
  int* indexArray,
  double* coefficientArray,
  long& position)
{
int Dim = particles->GetHilbertSpaceDimension();
double Coefficient;
int Index;

/* Diagonal part */
double TmpDiagonal = 0.0;

for (int m = 0; m <= this->LzMax; ++m)
  for (int a = 0; a < this->NbrFlavors; ++a)
    TmpDiagonal += this->OneBodyPotentials[a][a][m]
                   * particles->AdsigmaAsigma(index + this->PrecalculationShift, m, a);

if (TmpDiagonal != 0.0)
{
  indexArray[position] = index + this->PrecalculationShift;
  coefficientArray[position] = TmpDiagonal;
  ++position;
}

/* Off-diagonal */
for (int m = 0; m <= this->LzMax; ++m)
  for (int a = 0; a < this->NbrFlavors; ++a)
    for (int b = 0; b < this->NbrFlavors; ++b)
      if (a != b)
      {
        Index = particles->AdsigmaAsigma(index + this->PrecalculationShift,
                                         m, m, a, b, Coefficient);
        if (Index < Dim)
        {
          indexArray[position] = Index;
          coefficientArray[position] =
              Coefficient * this->OneBodyPotentials[a][b][m];
          ++position;
        }
      }
}

#endif
