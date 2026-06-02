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


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithNFlavor.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


/**************************************************************/
/* Default constructor                                       */
/**************************************************************/

FermionOnSphereWithNFlavor::FermionOnSphereWithNFlavor()
{
  this->HilbertSpaceDimension = 0;
  this->NbrFermions = 0;
  this->IncNbrFermions = 0;
  this->TotalLz = 0;
  this->LzMax = 0;
  this->NbrLzValue = 0;
  this->NbrFlavors = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;
  this->NbrFermionsPerFlavor = 0;
  this->MaximumLookUpShift = 0;
  this->LookUpTableMemorySize = 0;
  this->LookUpTableShift = 0;
  this->LookUpTable = 0;
  this->SignLookUpTable = 0;
  this->SignLookUpTableMask = 0;
  this->MaximumSignLookUp = 0;

  if (this->HilbertSpaceDimension <= 0)
{
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->SignLookUpTable = 0;
  this->SignLookUpTableMask = 0;
  return;
}
}


/**************************************************************/
/* Basic constructor                                         */
/**************************************************************/

FermionOnSphereWithNFlavor::FermionOnSphereWithNFlavor(int nbrFermions, int totalLz, int lzMax, int nbrFlavors, int* nbrPerFlavor, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = nbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = lzMax + 1;
  this->NbrFlavors = nbrFlavors;

  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();

  this->NbrFermionsPerFlavor = new int[nbrFlavors];

  int sum = 0;
  for (int i = 0; i < nbrFlavors; ++i)
  {
    this->NbrFermionsPerFlavor[i] = nbrPerFlavor[i];
    sum += nbrPerFlavor[i];
  }

  if (sum != nbrFermions)
  {
    cout << "Flavor populations do not sum to total particles!" << endl;
    exit(1);
  }

  long shiftedLz = (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1;

  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, shiftedLz);

  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;

  this->StateDescription = new unsigned long[this->HilbertSpaceDimension];
  this->StateHighestBit = new int[this->HilbertSpaceDimension];

  long tmpDim =
      this->GenerateStates(this->NbrFermions, this->LzMax, shiftedLz, 0l );

  if (tmpDim != this->HilbertSpaceDimension)
  {
    cout << "Mismatch in state counting!" << endl;
    exit(1);
  }

  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}


/**************************************************************/
/* Copy constructor (shared memory)                          */
/**************************************************************/

FermionOnSphereWithNFlavor::FermionOnSphereWithNFlavor(const FermionOnSphereWithNFlavor& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->NbrFlavors = fermions.NbrFlavors;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->NbrFermionsPerFlavor = fermions.NbrFermionsPerFlavor;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}


/**************************************************************/
/* Destructor                                                */
/**************************************************************/

FermionOnSphereWithNFlavor::~FermionOnSphereWithNFlavor()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
  {
    //int maxBit = this->StateHighestBit[0];
    int maxBit = (this->LzMax + 1) * this->NbrFlavors - 1;

    delete[] this->StateDescription;
    delete[] this->StateHighestBit;
    delete[] this->LookUpTableShift;

    for (int i = 0; i <= maxBit; ++i)
      delete[] this->LookUpTable[i];

    delete[] this->LookUpTable;
    delete[] this->SignLookUpTable;
    delete[] this->SignLookUpTableMask;
    delete[] this->NbrFermionsPerFlavor;
  }
}


/**************************************************************/
/* Assignment                                                */
/**************************************************************/

FermionOnSphereWithNFlavor&
FermionOnSphereWithNFlavor::operator=(const FermionOnSphereWithNFlavor& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->NbrFlavors = fermions.NbrFlavors;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->NbrFermionsPerFlavor = fermions.NbrFermionsPerFlavor;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension =
      (long) this->HilbertSpaceDimension;

  return *this;
}


/**************************************************************/
/* Clone                                                     */
/**************************************************************/

AbstractHilbertSpace* FermionOnSphereWithNFlavor::Clone()
{
  return new FermionOnSphereWithNFlavor(*this);
}


/**************************************************************/
/* Quantum number interface                                  */
/**************************************************************/

List<AbstractQuantumNumber*> FermionOnSphereWithNFlavor::GetQuantumNumbers()
{
  List<AbstractQuantumNumber*> L;
  L += new SzQuantumNumber(this->TotalLz);
  return L;
}

AbstractQuantumNumber* FermionOnSphereWithNFlavor::GetQuantumNumber(int)
{
  return new SzQuantumNumber(this->TotalLz);
}

AbstractHilbertSpace* FermionOnSphereWithNFlavor::ExtractSubspace(AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  return 0;
}


// ---- GENERIC OPERATORS ----
// Occupation operator a^+_m_sigma a_m_sigma
// Single creation-annihilation operator a^+_m_sigma1 a_n_sigma2
// Double annihilation operator a_m_sigma a_n_sigma

double FermionOnSphereWithNFlavor::AdsigmaAsigma(int index, int m, int sigma)
{
    int bit = m * this->NbrFlavors + sigma;

    if (this->StateDescription[index] & (0x1ul << bit))
        return 1.0;
    else
        return 0.0;
}

int FermionOnSphereWithNFlavor::AdsigmaAsigma(int index, int m, int n, int sigma_m, int sigma_n, double& coefficient)
{
    int create_bit = m * this->NbrFlavors + sigma_m;
    int annihilate_bit = n * this->NbrFlavors + sigma_n;
    // Guard: if annihilation bit not occupied, return immediately
    if ((this->StateDescription[index] & (0x1ul << annihilate_bit)) == 0x0ul)
    {
        coefficient = 0.0;
        return this->HilbertSpaceDimension;
    }
    return this->GenericAdA(index, create_bit, annihilate_bit, coefficient);
}



/*
double FermionOnSphereWithNFlavor::AsigmaAsigma(int index, int n1, int n2, int sigma1, int sigma2)
{
    this->ProdATemporaryState = this->StateDescription[index];

    int bit1 = n1 * this->NbrFlavors + sigma1;
    int bit2 = n2 * this->NbrFlavors + sigma2;

    if (((this->ProdATemporaryState & (0x1ul << bit1)) == 0) ||
        ((this->ProdATemporaryState & (0x1ul << bit2)) == 0) || (bit1 == bit2))
        return 0.0;

    this->ProdALzMax = this->StateHighestBit[index];

    // ===================  Annihilate at bit2  ===================
    double coeff = 1.0;
    for (int offset = 0; offset < 64; offset += 16)
    {
        int pos = bit2 + offset;
        if (pos >= 64) break;
        coeff *= this->SignLookUpTable[(this->ProdATemporaryState >> pos) & this->SignLookUpTableMask[pos]];
    }
    this->ProdATemporaryState &= ~(0x1ul << bit2);

    // ===================  Annihilate at bit1  ===================
    for (int offset = 0; offset < 64; offset += 16)
    {
        int pos = bit1 + offset;
        if (pos >= 64) break;
        coeff *= this->SignLookUpTable[(this->ProdATemporaryState >> pos) & this->SignLookUpTableMask[pos]];
    }
    this->ProdATemporaryState &= ~(0x1ul << bit1);

    if (this->ProdATemporaryState != 0x0ul)
    {
        while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
            --this->ProdALzMax;
    }
    else
        this->ProdALzMax = 0;

    return coeff;
}*/

double FermionOnSphereWithNFlavor::AsigmaAsigma(
  int index, int n1, int n2,
  int sigma1, int sigma2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  int bit1 = n1 * this->NbrFlavors + sigma1;
  int bit2 = n2 * this->NbrFlavors + sigma2;

  if (((this->ProdATemporaryState & (0x1ul << bit1)) == 0) ||
      ((this->ProdATemporaryState & (0x1ul << bit2)) == 0) ||
      (bit1 == bit2))
      return 0.0;

  this->ProdALzMax = this->StateHighestBit[index];

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> bit2) & this->SignLookUpTableMask[bit2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (bit2 + 16)) & this->SignLookUpTableMask[bit2 + 16]];

#ifdef __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (bit2 + 32)) & this->SignLookUpTableMask[bit2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (bit2 + 48)) & this->SignLookUpTableMask[bit2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << bit2);

  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> bit1) & this->SignLookUpTableMask[bit1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (bit1 + 16)) & this->SignLookUpTableMask[bit1 + 16]];

#ifdef __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (bit1 + 32)) & this->SignLookUpTableMask[bit1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (bit1 + 48)) & this->SignLookUpTableMask[bit1 + 48]];
#endif

  this->ProdATemporaryState &= ~(0x1ul << bit1);

  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
      --this->ProdALzMax;

  return Coefficient;
}


int FermionOnSphereWithNFlavor::AdsigmaAdsigma(int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
    unsigned long tmpState = this->ProdATemporaryState;
    int bit1 = m1 * this->NbrFlavors + sigma1;
    int bit2 = m2 * this->NbrFlavors + sigma2;

    if (((tmpState & (0x1ul << bit1)) != 0) ||
        ((tmpState & (0x1ul << bit2)) != 0) || (bit1 == bit2))
        return this->HilbertSpaceDimension;

    int newLzMax = this->ProdALzMax;

    // ===================  Create at bit2  =======================
    coefficient = 1.0;
    if (bit2 > newLzMax) newLzMax = bit2;
    /*else
    {
        for (int offset = 0; offset < 64; offset += 16)
        {
            int pos = bit2 + offset;
            if (pos >= 64) break;
            coefficient *= this->SignLookUpTable[(tmpState >> pos) & this->SignLookUpTableMask[pos]];
        }
    }*/
    else
    {
        coefficient *= this->SignLookUpTable[(tmpState >> bit2) & this->SignLookUpTableMask[bit2]];
        coefficient *= this->SignLookUpTable[(tmpState >> (bit2 + 16)) & this->SignLookUpTableMask[bit2 + 16]];
#ifdef __64_BITS__
        coefficient *= this->SignLookUpTable[(tmpState >> (bit2 + 32)) & this->SignLookUpTableMask[bit2 + 32]];
        coefficient *= this->SignLookUpTable[(tmpState >> (bit2 + 48)) & this->SignLookUpTableMask[bit2 + 48]];
#endif
    }
    tmpState |= (0x1ul << bit2);

    // ===================  Create at bit1  ======================
    if (bit1 > newLzMax) newLzMax = bit1;
    /*else
    {
        for (int offset = 0; offset < 64; offset += 16)
        {
            int pos = bit1 + offset;
            if (pos >= 64) break;
            coefficient *= this->SignLookUpTable[(tmpState >> pos) & this->SignLookUpTableMask[pos]];
        }
    }*/
    else
    {
        coefficient *= this->SignLookUpTable[(tmpState >> bit1) & this->SignLookUpTableMask[bit1]];
        coefficient *= this->SignLookUpTable[(tmpState >> (bit1 + 16)) & this->SignLookUpTableMask[bit1 + 16]];
#ifdef __64_BITS__
        coefficient *= this->SignLookUpTable[(tmpState >> (bit1 + 32)) & this->SignLookUpTableMask[bit1 + 32]];
        coefficient *= this->SignLookUpTable[(tmpState >> (bit1 + 48)) & this->SignLookUpTableMask[bit1 + 48]];
#endif
    }
    tmpState |= (0x1ul << bit1);
    return this->FindStateIndex(tmpState, newLzMax);
}




int FermionOnSphereWithNFlavor::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    {
      return this->HilbertSpaceDimension;
    }
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;

  // int TableMaxBit = (this->LzMax + 1) * this->NbrFlavors - 1;

  // if (lzmax < 0 || lzmax > TableMaxBit)
  //   return this->HilbertSpaceDimension;

  // unsigned long CurrentState = stateDescription >> this->LookUpTableShift[lzmax];

  // if (CurrentState > (unsigned long)this->LookUpTableMemorySize)
  //   CurrentState = this->LookUpTableMemorySize;

  // int PosMin = this->LookUpTable[lzmax][CurrentState];
  // int PosMax = this->LookUpTable[lzmax][CurrentState + 1];
  // int PosMid = (PosMin + PosMax) >> 1;

  // CurrentState = this->StateDescription[PosMid];

  // while ((PosMax != PosMid) && (CurrentState != stateDescription))
  // {
  //   if (CurrentState > stateDescription)
  //     PosMax = PosMid;
  //   else
  //     PosMin = PosMid;

  //   PosMid = (PosMin + PosMax) >> 1;
  //   CurrentState = this->StateDescription[PosMid];
  // }

  // if (CurrentState == stateDescription)
  //   return PosMid;
  // else
  //   return this->HilbertSpaceDimension;
}


ostream& FermionOnSphereWithNFlavor::PrintState(ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << " | ";

  for (int i = this->NbrLzValue - 1; i >= 0; --i)
  {
    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    {
      int bit = i * this->NbrFlavors + sigma;
      if (TmpState & (0x1ul << bit))
        Str << "f" << sigma << " ";
      else
        Str << "0 ";
    }
    Str << "| ";
  }

  return Str;
}
/*
long FermionOnSphereWithNFlavor::GenerateStates(int nbrFermions, int lzMax, int totalLz, long pos)
{
  if (nbrFermions < 0 || totalLz < 0)
    return pos;

  if (nbrFermions == 0)
  {
    if (totalLz == 0)
    {
      this->StateDescription[pos] = 0x0ul;
      return pos + 1;
    }
    return pos;
  }

  if (lzMax < 0)
    return pos;

  // Case 1: no particle in this orbital
  pos = this->GenerateStates(nbrFermions, lzMax - 1, totalLz, pos);

  // Case 2: choose subset of flavors occupying this orbital
  int maxMask = 1 << this->NbrFlavors;

  for (int mask = 1; mask < maxMask; ++mask)
  {
    int k = 0;
    bool allowed = true;

    // check flavor availability
    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    {
      if (mask & (1 << sigma))
      {
        if (this->NbrFermionsPerFlavor[sigma] <= 0)
        {
          allowed = false;
          break;
        }
        k++;
      }
    }

    if (!allowed)
      continue;

    if (k > nbrFermions)
      continue;

    if (totalLz - k * lzMax < 0)
      continue;

    // subtract flavor populations
    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
      if (mask & (1 << sigma))
        this->NbrFermionsPerFlavor[sigma]--;

    long tmpPos = this->GenerateStates(nbrFermions - k, lzMax - 1, totalLz - k * lzMax, pos);

    unsigned long bitmask = 0x0ul;

    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    {
      if (mask & (1 << sigma))
      {
        int bit = lzMax * this->NbrFlavors + sigma;
        bitmask |= (0x1ul << bit);
      }
    }

    for (; pos < tmpPos; ++pos)
      this->StateDescription[pos] |= bitmask;

    // restore flavor populations
    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
      if (mask & (1 << sigma))
        this->NbrFermionsPerFlavor[sigma]++;
  }

  return pos;
}*/

long FermionOnSphereWithNFlavor::GenerateStates(int nbrFermions, int lzMax, int totalLz, long pos)
{
  if (nbrFermions < 0 || totalLz < 0)
    return pos;

  if (nbrFermions == 0)
  {
    if (totalLz == 0)
    {
      this->StateDescription[pos] = 0x0ul;
      return pos + 1;
    }
    return pos;
  }

  if (lzMax < 0)
    return pos;

  int maxMask = 1 << this->NbrFlavors;

  // --------------------------------------------------
  // IMPORTANT:
  // iterate occupied masks first, in descending order,
  // and only then do the empty-orbital branch.
  // This matches the SU3/SU4 style ordering.
  // --------------------------------------------------
  for (int mask = maxMask - 1; mask >= 1; --mask)
  {
    int k = 0;
    bool allowed = true;

    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    {
      if (mask & (1 << sigma))
      {
        if (this->NbrFermionsPerFlavor[sigma] <= 0)
        {
          allowed = false;
          break;
        }
        ++k;
      }
    }

    if (!allowed)
      continue;

    if (k > nbrFermions)
      continue;

    if (totalLz - k * lzMax < 0)
      continue;

    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
      if (mask & (1 << sigma))
        --this->NbrFermionsPerFlavor[sigma];

    long tmpPos = this->GenerateStates(nbrFermions - k, lzMax - 1, totalLz - k * lzMax, pos);

    unsigned long bitmask = 0x0ul;
    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    {
      if (mask & (1 << sigma))
      {
        int bit = lzMax * this->NbrFlavors + sigma;
        bitmask |= (0x1ul << bit);
      }
    }

    for (; pos < tmpPos; ++pos)
      this->StateDescription[pos] |= bitmask;

    for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
      if (mask & (1 << sigma))
        ++this->NbrFermionsPerFlavor[sigma];
  }

  // empty orbital LAST, like SU3
  pos = this->GenerateStates(nbrFermions, lzMax - 1, totalLz, pos);

  return pos;
}


long FermionOnSphereWithNFlavor::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  int totalRemaining = 0;
  for (int s = 0; s < this->NbrFlavors; ++s)
    totalRemaining += this->NbrFermionsPerFlavor[s];

  if (nbrFermions < 0 || totalLz < 0)
    return 0l;
  
  if (nbrFermions == 0) return (totalLz == 0) ? 1l : 0l;
  /*{
    if (totalLz == 0)
      return 1l;
    return 0l;
  }*/

  if (lzMax < 0 || totalLz < 0)
    return 0l;

  if (totalRemaining != nbrFermions)
    return 0l;

  long Dim = 0l;
  int maxMask = 1 << this->NbrFlavors;

  // occupied masks first
  for (int mask = maxMask - 1; mask >= 1; --mask)
  {
    int k = 0;
    bool allowed = true;

    for (int s = 0; s < this->NbrFlavors; ++s)
    {
      if (mask & (1 << s))
      {
        if (this->NbrFermionsPerFlavor[s] <= 0)
        {
          allowed = false;
          break;
        }
        ++k;
      }
    }

    if (!allowed) continue;
    if (k > nbrFermions) continue;
    if (totalLz - k * lzMax < 0) continue;

    for (int s = 0; s < this->NbrFlavors; ++s)
      if (mask & (1 << s))
        --this->NbrFermionsPerFlavor[s];

    Dim += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - k, lzMax - 1, totalLz - k * lzMax);

    for (int s = 0; s < this->NbrFlavors; ++s)
      if (mask & (1 << s))
        ++this->NbrFermionsPerFlavor[s];
  }

  // empty orbital last
  Dim += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);
  
  /* // Case 1: no particle in this orbital
  Dim += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);

  // Case 2: choose subset of flavors to occupy this orbital
  int maxMask = 1 << this->NbrFlavors;
  for (int mask = 1; mask < maxMask; ++mask)
  {
    int k = 0;
    bool allowed = true;

    for (int s = 0; s < this->NbrFlavors; ++s)
    {
      if (mask & (1 << s))
      {
        if (this->NbrFermionsPerFlavor[s] <= 0)
        {
          allowed = false;
          break;
        }
        k++;
      }
    }

    if (!allowed)
      continue;

    if (totalLz - k * lzMax < 0)
      continue;

    // subtract
    for (int s = 0; s < this->NbrFlavors; ++s)
      if (mask & (1 << s))
        this->NbrFermionsPerFlavor[s]--;

    Dim += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz - k * lzMax);

    // restore
    for (int s = 0; s < this->NbrFlavors; ++s)
      if (mask & (1 << s))
        this->NbrFermionsPerFlavor[s]++;
  }*/

  return Dim;
}


void FermionOnSphereWithNFlavor::GenerateLookUpTable(unsigned long memory)
{
  // --------------------------------------------------------
  // 1. Compute highest set bit for each basis state
  // --------------------------------------------------------
  int MaxHighestBit = 0;

  if (this->HilbertSpaceDimension <= 0)
    return;

  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
  {
    unsigned long TmpPosition = this->StateDescription[i];

    int CurrentHighestBit = 0;

    if (TmpPosition != 0x0ul)
    {
#ifdef __64_BITS__
      CurrentHighestBit = 63;
#else
      CurrentHighestBit = 31;
#endif
      while ((CurrentHighestBit > 0) &&
             ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul))
        --CurrentHighestBit;
    }

    this->StateHighestBit[i] = CurrentHighestBit;

    if (CurrentHighestBit > MaxHighestBit)
      MaxHighestBit = CurrentHighestBit;
  }

  // full operator-accessible bit range
  int TheoreticalMaxBit = (this->LzMax + 1) * this->NbrFlavors - 1;
  int TableMaxBit = TheoreticalMaxBit;

  if (MaxHighestBit > TheoreticalMaxBit)
  {
    cout << "Warning: MaxHighestBit (" << MaxHighestBit
         << ") > TheoreticalMaxBit (" << TheoreticalMaxBit << ")"
         << endl;
  }

  // --------------------------------------------------------
  // 2. Evaluate lookup table size
  // --------------------------------------------------------
  memory /= (sizeof(int*) * (TableMaxBit + 1));

  this->MaximumLookUpShift = 1;
  while (memory > 0)
  {
    memory >>= 1;
    ++this->MaximumLookUpShift;
  }

  if (this->MaximumLookUpShift > TableMaxBit)
    this->MaximumLookUpShift = TableMaxBit;

  if (this->MaximumLookUpShift < 0)
    this->MaximumLookUpShift = 0;

  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // --------------------------------------------------------
  // 3. Allocate lookup tables for FULL theoretical range
  // --------------------------------------------------------
  this->LookUpTable = new int*[TableMaxBit + 1];
  this->LookUpTableShift = new int[TableMaxBit + 1];

  for (int i = 0; i <= TableMaxBit; ++i)
  {
    this->LookUpTable[i] = new int[this->LookUpTableMemorySize + 1];    
    this->LookUpTableShift[i] = 0;
    for (int j = 0; j <= this->LookUpTableMemorySize; ++j)
      this->LookUpTable[i][j] = 0;
  }

  // track which highest-bit sectors are actually present
  bool* BitPresent = new bool[TableMaxBit + 1];
  for (int i = 0; i <= TableMaxBit; ++i)
    BitPresent[i] = false;

  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
  {
    int hb = this->StateHighestBit[i];
    if ((hb >= 0) && (hb <= TableMaxBit))
      BitPresent[hb] = true;
  }

  // --------------------------------------------------------
  // 4. Fill lookup tables for sectors that actually appear
  // --------------------------------------------------------
  int CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];

  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] =
      CurrentHighestBit + 1 - this->MaximumLookUpShift;

  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue =
    this->StateDescription[0] >> CurrentShift;

  while (CurrentLookUpTableValue > TmpLookUpTableValue)
  {
    TmpLookUpTable[CurrentLookUpTableValue] = 0;
    --CurrentLookUpTableValue;
  }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;

  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
  {
    if (CurrentHighestBit != this->StateHighestBit[i])
    {
      while (CurrentLookUpTableValue > 0)
      {
        TmpLookUpTable[CurrentLookUpTableValue] = i;
        --CurrentLookUpTableValue;
      }
      TmpLookUpTable[0] = i;

      CurrentHighestBit = this->StateHighestBit[i];
      TmpLookUpTable = this->LookUpTable[CurrentHighestBit];

      if (CurrentHighestBit < this->MaximumLookUpShift)
        this->LookUpTableShift[CurrentHighestBit] = 0;
      else
        this->LookUpTableShift[CurrentHighestBit] =
          CurrentHighestBit + 1 - this->MaximumLookUpShift;

      CurrentShift = this->LookUpTableShift[CurrentHighestBit];
      TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
      CurrentLookUpTableValue = this->LookUpTableMemorySize;

      while (CurrentLookUpTableValue > TmpLookUpTableValue)
      {
        TmpLookUpTable[CurrentLookUpTableValue] = i;
        --CurrentLookUpTableValue;
      }
      TmpLookUpTable[CurrentLookUpTableValue] = i;
    }
    else
    {
      TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;

      if (TmpLookUpTableValue != CurrentLookUpTableValue)
      {
        while (CurrentLookUpTableValue > TmpLookUpTableValue)
        {
          TmpLookUpTable[CurrentLookUpTableValue] = i;
          --CurrentLookUpTableValue;
        }
        TmpLookUpTable[CurrentLookUpTableValue] = i;
      }
    }
  }

  while (CurrentLookUpTableValue > 0)
  {
    TmpLookUpTable[CurrentLookUpTableValue] =
      this->HilbertSpaceDimension - 1;
    --CurrentLookUpTableValue;
  }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;

  // --------------------------------------------------------
  // 5. Fill missing highest-bit sectors by copying nearest
  //    valid lower sector, or nearest higher one if needed
  // --------------------------------------------------------
  for (int b = 0; b <= TableMaxBit; ++b)
  {
    if (BitPresent[b])
      continue;

    int CopyFrom = -1;

    for (int k = b - 1; k >= 0; --k)
    {
      if (BitPresent[k])
      {
        CopyFrom = k;
        break;
      }
    }

    if (CopyFrom < 0)
    {
      for (int k = b + 1; k <= TableMaxBit; ++k)
      {
        if (BitPresent[k])
        {
          CopyFrom = k;
          break;
        }
      }
    }

    if (CopyFrom >= 0)
    {
      this->LookUpTableShift[b] = this->LookUpTableShift[CopyFrom];
      for (int j = 0; j <= this->LookUpTableMemorySize; ++j)
        this->LookUpTable[b][j] = this->LookUpTable[CopyFrom][j];
    }
  }

  delete[] BitPresent;

  // --------------------------------------------------------
  // 6. Sign lookup tables
  // --------------------------------------------------------
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double[Size];

  for (int j = 0; j < Size; ++j)
  {
    int Count = 0;
    int TmpNbr = j;
    while (TmpNbr != 0)
    {
      if (TmpNbr & 0x1) ++Count;
      TmpNbr >>= 1;
    }
    this->SignLookUpTable[j] = (Count & 1) ? -1.0 : 1.0;
  }

#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long[128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long[64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}

/*************************************************************
void FermionOnSphereWithNFlavor::GenerateLookUpTable(unsigned long memory)
{
  // determine highest bit for each state
  unsigned long TmpPosition = this->StateDescription[0];

#ifdef __64_BITS__
  int CurrentHighestBit = 63;
#else
  int CurrentHighestBit = 31;
#endif

  while (CurrentHighestBit > 0 && (TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;

  int MaxHighestBit = CurrentHighestBit;
  this->StateHighestBit[0] = CurrentHighestBit;

  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
  {
    TmpPosition = this->StateDescription[i];

    while (CurrentHighestBit > 0 && (TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
      --CurrentHighestBit;

    this->StateHighestBit[i] = CurrentHighestBit;
  }

  // evaluate look-up table size
  memory /= (sizeof(int*) * (MaxHighestBit + 1));

  this->MaximumLookUpShift = 1;
  while (memory > 0)
  {
    memory >>= 1;
    ++this->MaximumLookUpShift;
  }

  if (this->MaximumLookUpShift > MaxHighestBit)
    this->MaximumLookUpShift = MaxHighestBit;

  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // allocate lookup tables
  this->LookUpTable = new int*[MaxHighestBit + 1];
  this->LookUpTableShift = new int[MaxHighestBit + 1];

  for (int i = 0; i <= MaxHighestBit; ++i)
    this->LookUpTable[i] = new int[this->LookUpTableMemorySize + 1];

  CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];

  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;

  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;

  while (CurrentLookUpTableValue > TmpLookUpTableValue)
  {
    TmpLookUpTable[CurrentLookUpTableValue] = 0;
    --CurrentLookUpTableValue;
  }

  TmpLookUpTable[CurrentLookUpTableValue] = 0;

  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
  {
    if (CurrentHighestBit != this->StateHighestBit[i])
    {
      while (CurrentLookUpTableValue > 0)
      {
        TmpLookUpTable[CurrentLookUpTableValue] = i;
        --CurrentLookUpTableValue;
      }

      TmpLookUpTable[0] = i;

      CurrentHighestBit = this->StateHighestBit[i];
      TmpLookUpTable = this->LookUpTable[CurrentHighestBit];

      if (CurrentHighestBit < this->MaximumLookUpShift)
        this->LookUpTableShift[CurrentHighestBit] = 0;
      else
        this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;

      CurrentShift = this->LookUpTableShift[CurrentHighestBit];
      TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
      CurrentLookUpTableValue = this->LookUpTableMemorySize;

      while (CurrentLookUpTableValue > TmpLookUpTableValue)
      {
        TmpLookUpTable[CurrentLookUpTableValue] = i;
        --CurrentLookUpTableValue;
      }

      TmpLookUpTable[CurrentLookUpTableValue] = i;
    }
    else
    {
      TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;

      if (TmpLookUpTableValue != CurrentLookUpTableValue)
      {
        while (CurrentLookUpTableValue > TmpLookUpTableValue)
        {
          TmpLookUpTable[CurrentLookUpTableValue] = i;
          --CurrentLookUpTableValue;
        }

        TmpLookUpTable[CurrentLookUpTableValue] = i;
      }
    }
  }

  while (CurrentLookUpTableValue > 0)
  {
    TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
    --CurrentLookUpTableValue;
  }

  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;

  // sign look-up table
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double[Size];

  for (int j = 0; j < Size; ++j)
  {
    int Count = 0;
    int TmpNbr = j;

    while (TmpNbr != 0)
    {
      if (TmpNbr & 0x1)
        ++Count;
      TmpNbr >>= 1;
    }

    this->SignLookUpTable[j] = (Count & 1) ? -1.0 : 1.0;
  }

#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long[128];

  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;

  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);

  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long[64];

  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;

  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);

  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}*/


void FermionOnSphereWithNFlavor::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis)
{
  int* momentum = new int[this->NbrFermions];
  int* initialFlavor = new int[this->NbrFermions];
  int* currentFlavor = new int[this->NbrFermions];
  int* flavorCount = new int[this->NbrFlavors];

  targetState.ClearVector();

  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
  {
    unsigned long state = this->StateDescription[i];
    int count = 0;

    for (int s = 0; s < this->NbrFlavors; ++s)
      flavorCount[s] = 0;

    for (int m = this->LzMax; m >= 0; --m)
      for (int s = this->NbrFlavors - 1; s >= 0; --s)
      {
        int bit = m * this->NbrFlavors + s;
        if (state & (0x1ul << bit))
        {
          momentum[count] = m;
          initialFlavor[count] = s;
          count++;
        }
      }

    this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, momentum, initialFlavor, currentFlavor, flavorCount, oneBodyBasis);
  }

  delete[] momentum;
  delete[] initialFlavor;
  delete[] currentFlavor;
  delete[] flavorCount;
}


void FermionOnSphereWithNFlavor::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coeff, int pos, int* momentum, int* initialFlavor, int* currentFlavor, int* flavorCount, ComplexMatrix* oneBodyBasis)
{
  if (pos == this->NbrFermions)
  {
    for (int s = 0; s < this->NbrFlavors; ++s)
      if (flavorCount[s] != this->NbrFermionsPerFlavor[s])
        return;

    unsigned long newState = 0x0ul;
    unsigned long sign = 0x0ul;

    for (int i = 0; i < this->NbrFermions; ++i)
    {
      int bit = momentum[i] * this->NbrFlavors + currentFlavor[i];
      unsigned long mask = 0x1ul << bit;
      if (newState & mask)
        return;

      unsigned long tmp = newState & (mask - 1);
#ifdef __64_BITS__
      tmp ^= tmp >> 32;
#endif
      tmp ^= tmp >> 16; tmp ^= tmp >> 8; tmp ^= tmp >> 4; tmp ^= tmp >> 2;
      sign ^= (tmp ^ (tmp >> 1)) & 0x1ul;
      newState |= mask;
    }

    int highestBit = this->NbrLzValue * this->NbrFlavors - 1;
    while (highestBit > 0 && ((newState >> highestBit) == 0)) highestBit--;

    int idx = this->FindStateIndex(newState, highestBit);
    if (idx < this->HilbertSpaceDimension)
      targetState[idx] += (sign == 0 ? coeff : -coeff);

    return;
  }

  int m = momentum[pos];
  int oldS = initialFlavor[pos];

  for (int newS = 0; newS < this->NbrFlavors; ++newS)
  {
    if (flavorCount[newS] >= this->NbrFermionsPerFlavor[newS])
      continue;

    flavorCount[newS]++;
    currentFlavor[pos] = newS;

    this->TransformOneBodyBasisRecursive(targetState, coeff * oneBodyBasis[m][newS][oldS], pos + 1, momentum, initialFlavor, currentFlavor, flavorCount, oneBodyBasis);

    flavorCount[newS]--;
  }
}

