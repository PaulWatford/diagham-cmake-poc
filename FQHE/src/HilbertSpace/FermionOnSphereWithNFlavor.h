////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with N flavors               //
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


#ifndef FERMIONONSPHEREWITHNFLAVOR_H
#define FERMIONONSPHEREWITHNFLAVOR_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithNFlavor.h"

#include <iostream>
#include <cstring>
#include <cmath>

using std::cout;
using std::endl;

class FermionOnSphereWithNFlavor : public ParticleOnSphereWithNFlavor
{
 protected:

  int NbrFermions;
  int IncNbrFermions;

  int TotalLz;
  int LzMax;
  int NbrLzValue;

  //int NbrFlavors;
  int* NbrFermionsPerFlavor;

  int HighestBit;

  // array describing each state
  unsigned long* StateDescription;

  // maximum Lz bit for each state
  int* StateHighestBit;

  // lookup structures
  int MaximumLookUpShift;
  unsigned long LookUpTableMemorySize;
  int* LookUpTableShift;
  int** LookUpTable;

  double* SignLookUpTable;
  unsigned long* SignLookUpTableMask;
  int MaximumSignLookUp;

  // temporary state for A*A*
  unsigned long ProdATemporaryState;
  int ProdALzMax;

 public:

  FermionOnSphereWithNFlavor();

  virtual ~FermionOnSphereWithNFlavor();

  FermionOnSphereWithNFlavor(int nbrFermions, int totalLz, int lzMax, int nbrFlavors, int* nbrFermionsPerFlavor, unsigned long memory = 10000000);
  
  FermionOnSphereWithNFlavor(const FermionOnSphereWithNFlavor& fermions);
  
  FermionOnSphereWithNFlavor& operator=(const FermionOnSphereWithNFlavor& fermions);

  int GetParticleStatistic();

  virtual int GetNbrFlavors() const { return this->NbrFlavors; }

  virtual AbstractHilbertSpace* Clone();

  // Quantum number handling
virtual List<AbstractQuantumNumber*> GetQuantumNumbers();
virtual AbstractQuantumNumber* GetQuantumNumber(int);
virtual AbstractHilbertSpace* ExtractSubspace(AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

// Basis transformation
virtual void TransformOneBodyBasis(ComplexVector& initialState,
                                   ComplexVector& targetState,
                                   ComplexMatrix* oneBodyBasis);

virtual void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coeff, int pos, int* momentum,
                                            int* initialFlavor, int* currentFlavor, int* flavorCount,
                                            ComplexMatrix* oneBodyBasis);



  // ---- GENERIC OPERATORS ----

  double AdsigmaAsigma(int index, int m, int sigma);

  int AdsigmaAsigma(int index, int m, int n, int sigma_m, int sigma_n, double& coefficient);

  double AsigmaAsigma(int index, int n1, int n2, int sigma1, int sigma2);

  int AdsigmaAdsigma(int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // print
  ostream& PrintState (ostream& Str, int state);

 protected:

  // generalized bit mapping
  inline int BitIndex(int m, int sigma) const
  {
    return m * this->NbrFlavors + (this->NbrFlavors - 1 - sigma);
  }

  int GenericAdA(int index, int m, int n, double& coefficient);

  int FindStateIndex(unsigned long stateDescription, int lzmax);

  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  long GenerateStates(int nbrFermions, int lzMax, int totalLz, long pos);

  void GenerateLookUpTable(unsigned long memory);
};


inline int FermionOnSphereWithNFlavor::GetParticleStatistic()
{
  return AbstractQHEParticle::FermionicStatistic;
}


inline int FermionOnSphereWithNFlavor::GenericAdA(int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];

  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0x0ul))
  {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
  }

  int NewLargestBit = StateHighestBit;

  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];    // Is -ve sign needed?
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif

  State &= ~(0x1ul << n);

  if (State != 0x0ul)           //   if (NewLargestBit == n)
      while ((State >> NewLargestBit) == 0x0ul)
          --NewLargestBit;

  if ((State & (0x1ul << m)) != 0x0ul)
  {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
  }

  if (m > NewLargestBit)
      NewLargestBit = m;
  else
  {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
  }

  State |= (0x1ul << m);

  //return this->FindStateIndex(State, NewLargestBit);
  int newIndex = this->FindStateIndex(State, NewLargestBit);
  if (newIndex < 0 || newIndex >= this->HilbertSpaceDimension)
  {
    coefficient = 0.0;
    return this->HilbertSpaceDimension;
  }
  return newIndex;
}

/*
inline double FermionOnSphereWithNFlavor::AsigmaAsigma( int index, int n1, int n2, int sigma1, int sigma2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  int bit1 = n1 * this->NbrFlavors + (this->NbrFlavors - 1 - sigma1);
  int bit2 = n2 * this->NbrFlavors + (this->NbrFlavors - 1 - sigma2);

  if (((this->ProdATemporaryState & (0x1ul << bit1)) == 0) ||
      ((this->ProdATemporaryState & (0x1ul << bit2)) == 0) ||
      (bit1 == bit2))
      return 0.0;

  this->ProdALzMax = this->StateHighestBit[index];

  double Coefficient =
      this->SignLookUpTable[
          (this->ProdATemporaryState >> bit2)
          & this->SignLookUpTableMask[bit2]
      ];

  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> (bit2 + 16))
      & this->SignLookUpTableMask[bit2 + 16]
  ];

#ifdef __64_BITS__
  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> (bit2 + 32))
      & this->SignLookUpTableMask[bit2 + 32]
  ];
  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> (bit2 + 48))
      & this->SignLookUpTableMask[bit2 + 48]
  ];
#endif

  this->ProdATemporaryState &= ~(0x1ul << bit2);

  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> bit1)
      & this->SignLookUpTableMask[bit1]
  ];

  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> (bit1 + 16))
      & this->SignLookUpTableMask[bit1 + 16]
  ];

#ifdef __64_BITS__
  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> (bit1 + 32))
      & this->SignLookUpTableMask[bit1 + 32]
  ];
  Coefficient *= this->SignLookUpTable[
      (this->ProdATemporaryState >> (bit1 + 48))
      & this->SignLookUpTableMask[bit1 + 48]
  ];
#endif

  this->ProdATemporaryState &= ~(0x1ul << bit1);

  if (this->ProdATemporaryState != 0x0ul)
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
          --this->ProdALzMax;
  else
      this->ProdALzMax = 0;

  return Coefficient;
}


inline int FermionOnSphereWithNFlavor::AdsigmaAdsigma( int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
unsigned long TmpState = this->ProdATemporaryState;

// Convert (m, sigma) to global bit index
int bit1 = m1 * this->NbrFlavors + (this->NbrFlavors - 1 - sigma1);
int bit2 = m2 * this->NbrFlavors + (this->NbrFlavors - 1 - sigma2);

// If either state already occupied or same bit → forbidden
if (((TmpState & (0x1ul << bit1)) != 0) ||
  ((TmpState & (0x1ul << bit2)) != 0) ||
  (bit1 == bit2))
  return this->HilbertSpaceDimension;

int NewLzMax = this->ProdALzMax;

coefficient = 1.0;

// Insert bit2 first (SU4 ordering preserved)
if (bit2 > NewLzMax)
{
  NewLzMax = bit2;
}
else
{
  coefficient *= this->SignLookUpTable[
      (TmpState >> bit2) & this->SignLookUpTableMask[bit2]
  ];
  coefficient *= this->SignLookUpTable[
      (TmpState >> (bit2 + 16)) & this->SignLookUpTableMask[bit2 + 16]
  ];
#ifdef __64_BITS__
  coefficient *= this->SignLookUpTable[
      (TmpState >> (bit2 + 32)) & this->SignLookUpTableMask[bit2 + 32]
  ];
  coefficient *= this->SignLookUpTable[
      (TmpState >> (bit2 + 48)) & this->SignLookUpTableMask[bit2 + 48]
  ];
#endif
}

TmpState |= (0x1ul << bit2);

// Insert bit1
if (bit1 > NewLzMax)
{
  NewLzMax = bit1;
}
else
{
  coefficient *= this->SignLookUpTable[
      (TmpState >> bit1) & this->SignLookUpTableMask[bit1]
  ];
  coefficient *= this->SignLookUpTable[
      (TmpState >> (bit1 + 16)) & this->SignLookUpTableMask[bit1 + 16]
  ];
#ifdef __64_BITS__
  coefficient *= this->SignLookUpTable[
      (TmpState >> (bit1 + 32)) & this->SignLookUpTableMask[bit1 + 32]
  ];
  coefficient *= this->SignLookUpTable[
      (TmpState >> (bit1 + 48)) & this->SignLookUpTableMask[bit1 + 48]
  ];
#endif
}

TmpState |= (0x1ul << bit1);

return this->FindStateIndex(TmpState, NewLzMax);
}
*/


#endif
