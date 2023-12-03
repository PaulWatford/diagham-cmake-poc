////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                  class of particle on sphere with SU(12) spin              //
//                                                                            //
//                        last modification : 27/11/2023                      //
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
#include "HilbertSpace/FermionOnSphereWithSU12Spin.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

FermionOnSphereWithSU12Spin::FermionOnSphereWithSU12Spin ()
{
}

// destructor
//

FermionOnSphereWithSU12Spin::~FermionOnSphereWithSU12Spin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      unsigned long TmpPosition = this->StateDescription[0];
      int CurrentHighestBit = (this->LzMax + 1) * 12 - 1;
      while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i <= CurrentHighestBit; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSU12Spin::FindStateIndex(unsigned long stateDescription, int lzmax)
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
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereWithSU12Spin::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << " | ";
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i * 12)) & ((unsigned long) 0x3ful));
      for (int j = 0; j < 12; ++j)
	{
	  if (((Tmp >> j) & 0x1ul) != 0x0ul)
	    {
	      Str << (j + 1);
	    }
	  else
	    {
	      Str << " ";
	    }
	}
      Str << "| ";
    }
  return Str;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSU12Spin::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpPosition = this->StateDescription[0];
#ifdef __64_BITS__
  int CurrentHighestBit = 63;
#else
  int CurrentHighestBit = 31;
#endif
  while (((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul) && (CurrentHighestBit > 0))
    --CurrentHighestBit;
  
  if (this->StateHighestBit != 0)
    {
      this->StateHighestBit[0] = CurrentHighestBit;
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  TmpPosition = this->StateDescription[i];
	  while (((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul) && (CurrentHighestBit > 0))
	    --CurrentHighestBit;  
	  this->StateHighestBit[i] = CurrentHighestBit;
	}
      CurrentHighestBit = this->StateHighestBit[0];
    }

  // evaluate look-up table size
  memory /= (sizeof(int*) * 12 * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > (12 * this->NbrLzValue))
    this->MaximumLookUpShift = (12 * this->NbrLzValue);
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  CurrentHighestBit = this->StateHighestBit[0];
  this->LookUpTable = new int* [CurrentHighestBit + 1];
  this->LookUpTableShift = new int [CurrentHighestBit + 1];
  for (int i = 0; i <= CurrentHighestBit; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLargestBit = CurrentHighestBit;
  int* TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
  if (CurrentLargestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLargestBit] = 0;
  else
    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLargestBit];
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
      TmpPosition = this->StateDescription[i];
      while (((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul) && (CurrentHighestBit > 0))
	--CurrentHighestBit;  
      if (CurrentLargestBit != CurrentHighestBit)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  CurrentLargestBit--;
	  while (CurrentLargestBit > CurrentHighestBit)
	    {
	      if (CurrentLargestBit < this->MaximumLookUpShift)
		this->LookUpTableShift[CurrentLargestBit] = 0;
	      else
		this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
	      TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
	      CurrentLookUpTableValue = this->LookUpTableMemorySize;
	      while (CurrentLookUpTableValue > 0x0ul)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      CurrentLargestBit--;
	    }
	  TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
	  if (CurrentLargestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLargestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLargestBit];
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
  this->GenerateSignLookUpTable();
}

// generate look-up table for sign calculation
//

void FermionOnSphereWithSU12Spin::GenerateSignLookUpTable()
{
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}

