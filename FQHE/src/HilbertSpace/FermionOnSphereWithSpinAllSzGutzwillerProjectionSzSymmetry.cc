/////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                       all Sz sectors, Sz<->-Sz symmetry                    //
//                  and using Gutziller projection in orbital space           //
//                                                                            //
//                        last modification : 08/12/2023                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementMatrixOperation.h"
#include <math.h>
#include <bitset>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor 
//

FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// minusParity = select the Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry (int nbrFermions, int totalLz, int lzMax, bool minusParity, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
#ifdef __64_BITS__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 16 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 16 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 1;
    }
#endif

  // temporary space to accelerate state generation
  this->MaxTotalLz = new int*[2*NbrLzValue];
  for (int i=0; i<2*NbrLzValue; ++i)
    {
      MaxTotalLz[i] = new int[NbrFermions+1];
      for (int nbrFermions=0; nbrFermions<=NbrFermions; ++nbrFermions)
	{
	  MaxTotalLz[i][nbrFermions] = 0;
	  for (int n=0; (n<nbrFermions); ++n)	  
	    MaxTotalLz[i][nbrFermions] += (i-n)>>1;
 	}
    }

  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, (this->LzMax << 1) + 1 , (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, (this->LzMax << 1) + 1, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0l);


  this->SzParitySign = 1.0;
  if (minusParity == true)
    this->SzParitySign = -1.0;

  this->Flag.Initialize();
  this->TargetSpace = this;
  long TmpHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	{
	  this->StateDescription[i] = 0x0ul;
	}
      else
	{
	  this->GetStateSymmetry(this->StateDescription[i]);
	  if ((this->StateDescription[i] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
	    {
	      unsigned long TmpStateParity = this->StateDescription[i];
	      this->GetStateSingletParity(TmpStateParity);
	      if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (minusParity == false))
		  || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (minusParity == true)))
		{
		  ++TmpHilbertSpaceDimension;
		}
	      else
		{
		  this->StateDescription[i] = 0x0ul;
		}
	    }
	  else
	    {
	      ++TmpHilbertSpaceDimension;
	      this->StateDescription[i] &= ~FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT;
	    }
	}
    }
  cout << "dim = " << TmpHilbertSpaceDimension << endl;
  unsigned long* TmpStateDescription = new unsigned long [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0l;
  for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	  ++TmpHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) TmpHilbertSpaceDimension;

  if (this->LargeHilbertSpaceDimension > 0)
    {
      this->StateHighestBit =  new int [this->LargeHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      this->StateHighestBit = 0;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);

#ifdef __DEBUG__
      long UsedMemory = 0l;
      UsedMemory += this->LargeHilbertSpaceDimension * sizeof(unsigned long);
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024l)
	if (UsedMemory >= 1048576l)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
    }

#endif
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  // clean up temporary space
  for (int i=0; i<2*NbrLzValue; ++i)
    delete [] MaxTotalLz[i];
  delete [] MaxTotalLz;


}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->TargetSpace = this;
#ifdef __64_BITS__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 16 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 16 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 1;
    }
#endif
  if (this->HilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
      this->StateHighestBit = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
  UsedMemory = this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  cout << "memory requested for lookup table = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;

#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry(const FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LzParitySign = fermions.LzParitySign;
  this->SzParitySign = fermions.SzParitySign;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->TargetSpace = this;
}

// destructor
//

FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::~FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry& FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::operator = (const FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->LzParitySign = fermions.LzParitySign;
  this->SzParitySign = fermions.SzParitySign;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::Clone()
{
  return new FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry(*this);
}



// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::GenerateStates(int nbrFermions, int posMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (posMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = this->MaxTotalLz[posMax][nbrFermions];  
  
  if (LzTotalMax < totalLz)
    {
      return pos;
    }

  if (nbrFermions == 1)
    {
      if ((posMax>>1) >= totalLz)
	{
	  if ( ((posMax>>1) > totalLz) || (((posMax>>1) == totalLz) && (posMax&1)))
	    {
	      this->StateHighestBit[pos] = (totalLz << 1) + 1;
	      this->StateDescription[pos++] = 0x1ul << ((totalLz << 1) + 1);
	    }
	  this->StateHighestBit[pos] = (totalLz << 1);
	  this->StateDescription[pos++] = 0x1ul << (totalLz << 1);
	}
      return pos;
    }

  if (((posMax>>1) == 0) && (totalLz != 0))
    return pos;
  
  long TmpPos;
  unsigned long Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, posMax - 2, totalLz - (posMax>>1),  pos);
  Mask = 0x1ul << posMax;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, posMax - 2, totalLz - (posMax>>1),  pos);
  Mask = 0x1ul << (posMax - 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, posMax - 2, totalLz, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// posMax = highest position for next particle to be placed
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int posMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (posMax < (nbrFermions - 1)))
    return 0l;
  
  int LzTotalMax = this->MaxTotalLz[posMax][nbrFermions];

  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && ((posMax>>1) >= totalLz))
    {
      if ((posMax>>1) >totalLz)
	{
	  return 2l;
	}
      else
	{
	  return (1l + (posMax&1));
	}
    }
  long Count = 0l;
  Count += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, posMax - 2, totalLz - (posMax>>1));
  Count += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, posMax - 2, totalLz - (posMax>>1));
  Count += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, posMax - 2, totalLz);
  return Count;
}

// convert a given state from a generic basis to the gutzwiller basis
//
// state = reference on the vector to convert
// basis = pointer to the basis associated to state
// return value = converted vector

RealVector FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry::GutzwillerProjection(RealVector& state, ParticleOnSphereWithSpin* basis)
{
  FermionOnSphereWithSpinAllSzSzSymmetry* TmpSpace = (FermionOnSphereWithSpinAllSzSzSymmetry*) basis;
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NewLzMax = 1 + (this->LzMax << 1);
      unsigned long TmpState = this->StateDescription[i] & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
      int TmpIndex = TmpSpace->FindStateIndex(TmpState, NewLzMax);
      if (TmpIndex != TmpSpace->GetHilbertSpaceDimension())
	{
	  TmpVector[i] = state[TmpIndex];
	}
    }
  return TmpVector;
}
