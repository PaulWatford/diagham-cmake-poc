/////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                        all Sz sectors and Sz<->-Sz symmetry                //
//                                                                            //
//                        last modification : 12/07/2023                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinAllSzSzSymmetry.h"
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

FermionOnSphereWithSpinAllSzSzSymmetry::FermionOnSphereWithSpinAllSzSzSymmetry ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// minusParity = select the Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSzSzSymmetry::FermionOnSphereWithSpinAllSzSzSymmetry (int nbrFermions, int totalLz, int lzMax, bool minusParity, unsigned long memory)
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

FermionOnSphereWithSpinAllSzSzSymmetry::FermionOnSphereWithSpinAllSzSzSymmetry (char* fileName, unsigned long memory)
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

FermionOnSphereWithSpinAllSzSzSymmetry::FermionOnSphereWithSpinAllSzSzSymmetry(const FermionOnSphereWithSpinAllSzSzSymmetry& fermions)
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

FermionOnSphereWithSpinAllSzSzSymmetry::~FermionOnSphereWithSpinAllSzSzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinAllSzSzSymmetry& FermionOnSphereWithSpinAllSzSzSymmetry::operator = (const FermionOnSphereWithSpinAllSzSzSymmetry& fermions)
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

AbstractHilbertSpace* FermionOnSphereWithSpinAllSzSzSymmetry::Clone()
{
  return new FermionOnSphereWithSpinAllSzSzSymmetry(*this);
}


// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector FermionOnSphereWithSpinAllSzSzSymmetry::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpinAllSz& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (long i = 0l; i < nbodyBasis.GetLargeHilbertSpaceDimension(); ++i)
    {
      unsigned long TmpState = nbodyBasis.StateDescription[i];
      double Coefficient = 1.0;
      this->ProdASignature = 0x0ul;
      int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
      if (TmpPos != this->HilbertSpaceDimension)
	{
	  TmpVector[i] = state[TmpPos] * Coefficient;
	}
    }
  return TmpVector;  
}

// convert a given state from the usual n-body basis to the symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinAllSzSzSymmetry::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpinAllSz& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      if ((TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK) == Signature)
	{
	  Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == 0x0ul)
	    --NewLzMax;
	  if (Signature != 0x0ul)	
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += state[i] * M_SQRT1_2;
	  else
	    {
	      Signature = TmpState;
	      this->GetStateSingletParity(Signature);
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (this->LzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (this->LzParitySign < 0.0)))
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
	    }
	}
      else
	{
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == 0x0ul)
	    --NewLzMax;
	  Signature = TmpState;
	  this->GetStateSingletParity(Signature);
	  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->LzParitySign * state[i] * M_SQRT1_2;
	}
    }
  return TmpVector;  
}



// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSpinAllSzSzSymmetry::AduAd (int index, int m1, int n2, double& Coefficient)
{

  this->ProdATemporaryState = this->StateDescription[index];
  n2 <<= 1;
  unsigned long TmpMask = (0x1ul << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return this->HilbertSpaceDimension;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);

  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;

  if ((TmpState & (0x1ul << m1)) != 0x0ul) 
    return this->HilbertSpaceDimension;

  Coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  
  return this->SymmetrizeAdAdResult(TmpState, Coefficient);

}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinAllSzSzSymmetry::AddAu (int index, int m1, int n2, double& Coefficient)
{

  this->ProdATemporaryState = this->StateDescription[index];
  n2 <<= 1;
  ++n2;

  unsigned long TmpMask = (0x1ul << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return this->HilbertSpaceDimension;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);

  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;

  if ((TmpState & (0x1ul << m1)) != 0x0ul) 
    return this->HilbertSpaceDimension;

  Coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  
  return this->SymmetrizeAdAdResult(TmpState, Coefficient);
}
  
// Project the state from the tunneling space (all Sz's)
// to the space with the fixed projection of Sz (given by SzValue)
//
// state = state that needs to be projected
// su2Space = the subspace onto which the projection is carried out
// SzValue = the desired value of Sz

RealVector FermionOnSphereWithSpinAllSzSzSymmetry::ForgeSU2FromTunneling(RealVector& state, FermionOnSphereWithSpinSzSymmetry& su2Space, int SzValue)
{
  RealVector FinalState(su2Space.GetHilbertSpaceDimension(), true);
  int counter=0;
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      TmpState= this->GetSignedCanonicalState(TmpState);
	 TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      int TmpPos = this->LzMax << 1;
      int TmpSzValue=0;
      while (TmpPos >= 0)
	{
	  TmpSzValue+=((TmpState>>(TmpPos+1)) & 0x1ul) - ( (TmpState>>TmpPos) & 0x1ul );
	  TmpPos -= 2;
	}
      if (TmpSzValue == SzValue)
	{ 
	  ++counter;
	  int NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == 0x0ul)
	    --NewLzMax;
	  FinalState[su2Space.FindStateIndex(TmpState, NewLzMax)] += state[j];
	  cout<<"su2 "<< su2Space.FindStateIndex(TmpState, NewLzMax)<<" "<<TmpState<<endl;
	}
    }
  cout << "Nbr of stored components = " << counter << endl;
  //FinalState /= FinalState.Norm();
  return FinalState;  
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinAllSzSzSymmetry::GenerateStates(int nbrFermions, int posMax, int totalLz, long pos)
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
  TmpPos = this->GenerateStates(nbrFermions - 1, posMax - 1, totalLz - (posMax>>1),  pos);
  Mask = 0x1ul << posMax;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, posMax - 1, totalLz, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// posMax = highest position for next particle to be placed
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereWithSpinAllSzSzSymmetry::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int posMax, int totalLz)
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
  return  (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, posMax - 1, totalLz - (posMax>>1))
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, posMax - 1, totalLz));
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix FermionOnSphereWithSpinAllSzSzSymmetry::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0))
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	      TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector - (NbrFermionsComplementarySector * subsytemSize);

  long TmpNbrNonZeroElements = 0;

  if (nbrFermionSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0.0;
	  double Coefficient = 0.0;
 	  FermionOnSphereWithSpinAllSz TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize);
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
	      Coefficient = 1.0;
	      this->ProdASignature = 0x0ul;
	      int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpNbrNonZeroElements++;
		  TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex, Coefficient * groundState[TmpPos]);	
		}
            }	  

          if (TmpNbrNonZeroElements == 0)
            {
              RealMatrix TmpEntanglementMatrix;
              return TmpEntanglementMatrix;
            }

	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }


  if (NbrFermionsComplementarySector == 0)
    {
      FermionOnSphereWithSpinAllSz TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1, true);
      int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	    TmpEntanglementMatrix.AddToMatrixElement(i, 0, groundState[MinIndex + i]);
	}
      return TmpEntanglementMatrix;
    }

  FermionOnSphereWithSpinAllSz TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
 
  FermionOnSphereWithSpinAllSz TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize);
 
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  
  TmpNbrNonZeroElements = 0;
  double Coefficient = 0.0;

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	  Coefficient = 1.0;
	  this->ProdASignature = 0x0ul;
	  int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
              TmpNbrNonZeroElements++;
              TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, Coefficient * groundState[TmpPos]);
	    }
	}

     }

  if (TmpNbrNonZeroElements == 0)
   {
     RealMatrix TmpEntanglementMatrix;
     return TmpEntanglementMatrix;
   }
  return TmpEntanglementMatrix;    
}



// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsystem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSymmetrySector = Sz<->-Sz symmetry sector for particles that belong to the subsystem 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix FermionOnSphereWithSpinAllSzSzSymmetry::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, int szSymmetrySector, RealVector& groundState)
{
  int TotalSzSymmetrySector = 1;
  if (this->SzParitySign < 0.0)
    {
      TotalSzSymmetrySector = -1;
    }
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0) && (TotalSzSymmetrySector == 1.0))
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions) && (TotalSzSymmetrySector == szSymmetrySector))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	      TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector - (NbrFermionsComplementarySector * subsytemSize);

  long TmpNbrNonZeroElements = 0;

  if (nbrFermionSector == 0)
    {
      if ((lzSector == 0) && (szSymmetrySector == 1))
	{
	  double TmpValue = 0.0;
	  double Coefficient = 0.0;
 	  FermionOnSphereWithSpinAllSzSzSymmetry TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, (TotalSzSymmetrySector == -1));
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] & FERMION_SPHERE_SU2_SYMMETRIC_MASK) << (subsytemSize << 1);
	      Coefficient = 1.0;
	      this->ProdASignature = TmpHilbertSpace.StateDescription[MinIndex] & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
	      int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpNbrNonZeroElements++;
		  TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex,  groundState[TmpPos] * Coefficient);	
		}
            }	  

          if (TmpNbrNonZeroElements == 0)
            {
              RealMatrix TmpEntanglementMatrix;
              return TmpEntanglementMatrix;
            }

	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }


  if (NbrFermionsComplementarySector == 0)
    {
      if (szSymmetrySector == TotalSzSymmetrySector)
	{
	  FermionOnSphereWithSpinAllSzSzSymmetry TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1, (szSymmetrySector == -1));
	  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
	  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1, true);
	  int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
	  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(i, 0, groundState[MinIndex + i]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  FermionOnSphereWithSpinAllSzSzSymmetry TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1, (szSymmetrySector == -1));
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
 
  FermionOnSphereWithSpinAllSzSzSymmetry TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, ((szSymmetrySector * TotalSzSymmetrySector) == -1));
 
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  
  TmpNbrNonZeroElements = 0;
  double Coefficient = 0.0;

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = (TmpHilbertSpace.StateDescription[MinIndex] & FERMION_SPHERE_SU2_SYMMETRIC_MASK) << (subsytemSize << 1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = (TmpDestinationHilbertSpace.StateDescription[j] & FERMION_SPHERE_SU2_SYMMETRIC_MASK) | TmpComplementaryState;
	  if ((TmpDestinationHilbertSpace.StateDescription[j] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) != 0x0ul)
	    {
	      Coefficient = 1.0 / M_SQRT2;
	    }
	  else
	    {
	      Coefficient = M_SQRT2;//1.0;
	    }
	  if ((TmpHilbertSpace.StateDescription[MinIndex] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
	    {
	       Coefficient *= 1.0 / M_SQRT2;
	    }
	  this->ProdASignature = (TmpDestinationHilbertSpace.StateDescription[j] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT);//0x0ul;
	  int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
              TmpNbrNonZeroElements++;
              TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, Coefficient * groundState[TmpPos]);
	    }
	  if ((TmpDestinationHilbertSpace.StateDescription[j] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) != 0x0ul)
	    {
	      Coefficient = ((double) szSymmetrySector) / M_SQRT2;
	      if ((TmpHilbertSpace.StateDescription[MinIndex] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
		{
		  Coefficient *= 1.0 / M_SQRT2;
		}
	      TmpState = this->ApplySzSymmetry(TmpDestinationHilbertSpace.StateDescription[j] & FERMION_SPHERE_SU2_SYMMETRIC_MASK, Coefficient) | TmpComplementaryState;
	      this->ProdASignature = (TmpDestinationHilbertSpace.StateDescription[j] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT);//0x0ul;
	      int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpNbrNonZeroElements++;
		  TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, Coefficient * groundState[TmpPos]);
		}
	    }
	}

     }

  if (TmpNbrNonZeroElements == 0)
   {
     RealMatrix TmpEntanglementMatrix;
     return TmpEntanglementMatrix;
   }
  return TmpEntanglementMatrix;    
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix FermionOnSphereWithSpinAllSzSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, 
												       bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  int nbrOrbitalA = this->LzMax + 1;
  int nbrOrbitalB = this->LzMax + 1;  

  if (nbrParticleSector == 0)
    {
      if (lzSector == 0)
        {
          FermionOnSphereWithSpinAllSz TmpHilbertSpace(this->NbrFermions, this->TotalLz - lzSector, this->LzMax);
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
 	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
	      double Coefficient = 1.0;
	      this->ProdASignature = 0x0ul;
	      int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex, Coefficient * groundState[TmpPos]);	
		}
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;
        }
    }
 
 
  if (nbrParticleSector == this->NbrFermions)
    {
      if (lzSector == this->TotalLz)
        {
          FermionOnSphereWithSpinAllSz TmpDestinationHilbertSpace(nbrParticleSector, lzSector, this->LzMax);
          RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1,true);
	  for (int MinIndex = 0; MinIndex < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
            {
 	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[MinIndex];
	      double Coefficient = 1.0;
	      this->ProdASignature = 0x0ul;
	      int TmpPos = this->SymmetrizeAdAdResult(TmpState, Coefficient);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpEntanglementMatrix.AddToMatrixElement(MinIndex, 0, Coefficient * groundState[TmpPos]);	
		}
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;  
        }
    }
 
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryLzSector = this->TotalLz - lzSector;

  FermionOnSphereWithSpinAllSz SubsytemSpace(nbrParticleSector, lzSector, this->LzMax);
  FermionOnSphereWithSpinAllSz ComplementarySubsytemSpace(ComplementaryNbrParticles, this->TotalLz - lzSector, this->LzMax);

  if ((SubsytemSpace.GetHilbertSpaceDimension() > 0) && (ComplementarySubsytemSpace.GetHilbertSpaceDimension() > 0))
    {
      RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);
      
      long TmpNbrNonZeroElements = 0l;
      if (architecture != 0)
	{
	  FQHESphereParticleEntanglementMatrixOperation TmpOperation (this, &SubsytemSpace, &ComplementarySubsytemSpace, groundState, TmpEntanglementMatrix, removeBinomialCoefficient);
	  TmpOperation.ApplyOperation(architecture);
	  TmpNbrNonZeroElements = TmpOperation.GetNbrNonZeroMatrixElements();
	}
      else
	{
	  TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore(0, ComplementarySubsytemSpace.GetHilbertSpaceDimension(),
												&ComplementarySubsytemSpace, &SubsytemSpace, 
												groundState, &TmpEntanglementMatrix, removeBinomialCoefficient);
	}
      if (TmpNbrNonZeroElements > 0l)
	{
	  return TmpEntanglementMatrix;
	}
    }
  RealMatrix TmpEntanglementMatrixZero;
  return TmpEntanglementMatrixZero;
}
   
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz,Sz=0, Sz parity sectors.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated. It should be equal to zero
// szParitySector = parity sector for the discrete symmetry Sz<->-Sz
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix FermionOnSphereWithSpinAllSzSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, int szParity, RealVector& groundState, 
												       bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  int TotalSzSymmetrySector = 1;
  if (this->SzParitySign < 0.0)
    {
      TotalSzSymmetrySector = -1;
    }

  if (nbrParticleSector == 0)
    {
      if ((lzSector == 0) && (szParity== 1))
        {
          FermionOnSphereWithSpinAllSzSzSymmetry TmpHilbertSpace(this->NbrFermions, this->TotalLz - lzSector, this->LzMax, (TotalSzSymmetrySector == -1));
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              int TmpLzMax = (this->LzMax <<1) + 1;
              unsigned long TmpState = this->StateDescription[i];
              while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;             
              TmpEntanglementMatrix.SetMatrixElement(0, TmpHilbertSpace.FindStateIndex(TmpState, TmpLzMax), groundState[i]);
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;
        }
    }
 
 
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((lzSector == this->TotalLz) && (TotalSzSymmetrySector == szParity))
        {
          FermionOnSphereWithSpinAllSzSzSymmetry TmpDestinationHilbertSpace(nbrParticleSector, lzSector, this->LzMax, (TotalSzSymmetrySector == -1));
          RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1,true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              int TmpLzMax = (this->LzMax << 1) + 1;
              unsigned long TmpState = this->StateDescription[i];
              while ((TmpState >> TmpLzMax) == 0x0ul)
                --TmpLzMax;
	      
              TmpEntanglementMatrix.SetMatrixElement(TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax), 0, groundState[i]);
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;  
        }
    }
 
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryLzSector = this->TotalLz - lzSector;

  FermionOnSphereWithSpinAllSzSzSymmetry SubsytemSpace(nbrParticleSector, lzSector, this->LzMax, (szParity == -1));
  FermionOnSphereWithSpinAllSzSzSymmetry ComplementarySubsytemSpace(ComplementaryNbrParticles, this->TotalLz - lzSector, this->LzMax, ((szParity * TotalSzSymmetrySector) == -1));

  if ((SubsytemSpace.GetHilbertSpaceDimension() > 0) && (ComplementarySubsytemSpace.GetHilbertSpaceDimension() > 0))
    {
      RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);
      
      long TmpNbrNonZeroElements = 0l;
      if (architecture != 0)
	{
	  FQHESphereParticleEntanglementMatrixOperation TmpOperation (this, &SubsytemSpace, &ComplementarySubsytemSpace, groundState, TmpEntanglementMatrix, removeBinomialCoefficient);
	  TmpOperation.ApplyOperation(architecture);
	  TmpNbrNonZeroElements = TmpOperation.GetNbrNonZeroMatrixElements();
	}
      else
	{
	  TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore(0, ComplementarySubsytemSpace.GetHilbertSpaceDimension(),
												&ComplementarySubsytemSpace, &SubsytemSpace, 
												groundState, &TmpEntanglementMatrix, removeBinomialCoefficient);
	}
      if (TmpNbrNonZeroElements > 0l)
	{
	  return TmpEntanglementMatrix;
	}
    }
  RealMatrix TmpEntanglementMatrixZero;
  return TmpEntanglementMatrixZero;
}
   
// core part of the entanglement matrix evaluation for the particle partition
// 
// minIndex = first index to consider in the complementary Hilbert space
// nbrIndex = number of indices to consider in the complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// entanglementMatrix = pointer to entanglement matrix
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = number of components that have been added to the entanglement matrix

long FermionOnSphereWithSpinAllSzSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, 
												     ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealMatrix* entanglementMatrix, 
												     bool removeBinomialCoefficient)
{
  long TmpNbrNonZeroElements = 0l;
  FermionOnSphereWithSpinAllSz* TmpSubsystemSpace = (FermionOnSphereWithSpinAllSz*) destinationHilbertSpace;
  FermionOnSphereWithSpinAllSz* TmpComplementarySubsystemSpace = (FermionOnSphereWithSpinAllSz*) complementaryHilbertSpace;
  bool TmpSzParityFlag = true;
  FermionOnSphereWithSpinAllSz* SubsystemSpaceFull = new FermionOnSphereWithSpinAllSz(TmpSubsystemSpace->NbrFermions, TmpSubsystemSpace->TotalLz, TmpSubsystemSpace->LzMax);
  FermionOnSphereWithSpinAllSz* ComplementarySubsystemSpaceFull = new FermionOnSphereWithSpinAllSz(TmpComplementarySubsystemSpace->NbrFermions, TmpComplementarySubsystemSpace->TotalLz, TmpComplementarySubsystemSpace->LzMax);
  if ((SubsystemSpaceFull->GetHilbertSpaceDimension() == TmpSubsystemSpace->GetHilbertSpaceDimension()) &&
      (TmpComplementarySubsystemSpace->GetHilbertSpaceDimension() == ComplementarySubsystemSpaceFull->GetHilbertSpaceDimension()))
    {
      TmpSzParityFlag = false;
      delete SubsystemSpaceFull;
    }
  delete ComplementarySubsystemSpaceFull;
  if (TmpSzParityFlag == false)
    {
      FermionOnSphereWithSpinAllSz* SubsystemSpace = (FermionOnSphereWithSpinAllSz*) ((FermionOnSphereWithSpinAllSz*) destinationHilbertSpace)->Clone();
      FermionOnSphereWithSpinAllSz* ComplementarySubsystemSpace = (FermionOnSphereWithSpinAllSz*) ((FermionOnSphereWithSpinAllSz*) complementaryHilbertSpace)->Clone();
      
       double TmpInvBinomial = 1.0;
      if (removeBinomialCoefficient == false)
	{
	  BinomialCoefficients TmpBinomial (this->NbrFermions);
	  TmpInvBinomial = sqrt(1.0 / (TmpBinomial( this->NbrFermions , SubsystemSpace->NbrFermions))); 
	}

      int MaxIndex = minIndex + nbrIndex;
      double Coefficient = 0.0;
      for (int MinIndex = minIndex; MinIndex < MaxIndex; ++MinIndex)    
	{
	  unsigned long TmpState = ComplementarySubsystemSpace->StateDescription[MinIndex] & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  for (int j = 0; j < SubsystemSpace->HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = SubsystemSpace->StateDescription[j] & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	      if ((TmpState & TmpState2) == 0x0ul)
		{
		  unsigned long TmpState3 = TmpState | TmpState2;
		  this->ProdASignature = 0x0ul;
		  Coefficient = 1.0;
		  int TmpPos = this->SymmetrizeAdAdResult(TmpState3, Coefficient);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      Coefficient *= TmpInvBinomial;
		      unsigned long Sign = 0x0ul;
		      int Pos2 = (SubsystemSpace->LzMax << 1) + 1;
		      while ((Pos2 > 0) && (TmpState2 != 0x0ul))
			{
			  while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			    --Pos2;
			  TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
			  TmpState3 ^= TmpState3 >> 32;
#endif 
			  TmpState3 ^= TmpState3 >> 16;
			  TmpState3 ^= TmpState3 >> 8;
			  TmpState3 ^= TmpState3 >> 4;
			  TmpState3 ^= TmpState3 >> 2;
			  TmpState3 ^= TmpState3 >> 1;
			  Sign ^= TmpState3;
			  TmpState2 &= ~(0x1ul << Pos2);
			  --Pos2;
			}
		      if ((Sign & 0x1ul) == 0x0ul)           
			Coefficient *= 1.0;
		      else
			Coefficient *= -1.0;                 
		      TmpNbrNonZeroElements++;
		      entanglementMatrix->SetMatrixElement(j, MinIndex, Coefficient * groundState[TmpPos]);
		    }
		}
	    }
	}
      
      delete SubsystemSpace;
      delete ComplementarySubsystemSpace;
    }
  else
    {
      FermionOnSphereWithSpinAllSzSzSymmetry* SubsystemSpace = (FermionOnSphereWithSpinAllSzSzSymmetry*) ((FermionOnSphereWithSpinAllSzSzSymmetry*) destinationHilbertSpace)->Clone();
      FermionOnSphereWithSpinAllSzSzSymmetry* ComplementarySubsystemSpace = (FermionOnSphereWithSpinAllSzSzSymmetry*) ((FermionOnSphereWithSpinAllSzSzSymmetry*) complementaryHilbertSpace)->Clone();

      double TmpInvBinomial = 1.0;
      if (removeBinomialCoefficient == false)
	{
	  BinomialCoefficients TmpBinomial (this->NbrFermions);
	  TmpInvBinomial = sqrt(1.0 / (TmpBinomial( this->NbrFermions , SubsystemSpace->NbrFermions))); 
	}

      int MaxIndex = minIndex + nbrIndex;
      double Coefficient = 0.0;
      double Coefficient2 = 0.0;
      for (int MinIndex = minIndex; MinIndex < MaxIndex; ++MinIndex)    
	{
	  unsigned long TmpState = ComplementarySubsystemSpace->StateDescription[MinIndex] & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  double TmpRescalingFactor = 1.0;
	  if ((ComplementarySubsystemSpace->StateDescription[MinIndex] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) != 0x0ul)
	    {
	      TmpRescalingFactor = M_SQRT2;
	    }
	  for (int j = 0; j < SubsystemSpaceFull->HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = SubsystemSpaceFull->StateDescription[j] & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	      if ((TmpState & TmpState2) == 0x0ul)
		{
		  Coefficient2 = TmpRescalingFactor;
		  SubsystemSpace->ProdASignature = 0x0ul;
		  unsigned long TmpState4 = TmpState2;
		  int TmpPos2 =  SubsystemSpace->SymmetrizeAdAdResult(TmpState4, Coefficient2);
		  if (TmpPos2 != SubsystemSpace->HilbertSpaceDimension)
		    {
		      // cout << SubsystemSpace->NbrFermions << " " << SubsystemSpace->TotalLz  << endl;
		      // if ((SubsystemSpace->NbrFermions == 2) && (SubsystemSpace->TotalLz == -10))
		      // 	cout << MinIndex << " " << j << " " << Coefficient2 << " " << TmpPos2 << " " << SubsystemSpace->HilbertSpaceDimension << endl;
		      unsigned long TmpState3 = TmpState | TmpState2;
		      this->ProdASignature = 0x0ul;
		      Coefficient = 1.0;
		      int TmpPos = this->SymmetrizeAdAdResult(TmpState3, Coefficient);
		      Coefficient = Coefficient;
		      if (TmpPos != this->HilbertSpaceDimension)
			{
			  // if ((SubsystemSpace->NbrFermions == 2) && (SubsystemSpace->TotalLz == -10))
			  //   cout << Coefficient << endl;
			  Coefficient *= TmpInvBinomial;
			  unsigned long Sign = 0x0ul;
			  int Pos2 = (SubsystemSpace->LzMax << 1) + 1;
			  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
			    {
			      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
				--Pos2;
			      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
			      TmpState3 ^= TmpState3 >> 32;
#endif 
			      TmpState3 ^= TmpState3 >> 16;
			      TmpState3 ^= TmpState3 >> 8;
			      TmpState3 ^= TmpState3 >> 4;
			      TmpState3 ^= TmpState3 >> 2;
			      TmpState3 ^= TmpState3 >> 1;
			      Sign ^= TmpState3;
			      TmpState2 &= ~(0x1ul << Pos2);
			      --Pos2;
			    }
			  if ((Sign & 0x1ul) == 0x0ul)           
			    Coefficient *= 1.0;
			  else
			    Coefficient *= -1.0;                 
			  TmpNbrNonZeroElements++;
			  entanglementMatrix->AddToMatrixElement(TmpPos2, MinIndex,  groundState[TmpPos] * Coefficient * Coefficient2);
			}
		    }
		}
	    }
	}
      delete SubsystemSpaceFull;
      delete SubsystemSpace;
      delete ComplementarySubsystemSpace;
    }
  return TmpNbrNonZeroElements;
}
   

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsystem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& FermionOnSphereWithSpinAllSzSzSymmetry::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
																       int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
																       int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown,
																       RealMatrix& entanglementMatrix)
{
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrFermions, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrParticleSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      return entanglementMatrix;	  
    }
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
  
  FermionOnSphereWithSpinAllSz* SubsystemSpace = 0;
  int SubsystemSpaceDimension = 1;
  if (nbrParticleSector > 0)
    {      
      SubsystemSpace = new FermionOnSphereWithSpinAllSz(nbrParticleSector, lzSector, nbrOrbitalA - 1);
      SubsystemSpaceDimension = SubsystemSpace->GetHilbertSpaceDimension();
    }
  FermionOnSphereWithSpinAllSz* ComplementarySubsystemSpace = 0;
  int ComplementarySubsystemSpaceDimension = 1;  
  if (ComplementaryNbrParticles > 0)
    {
      ComplementarySubsystemSpace = new FermionOnSphereWithSpinAllSz(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1);
      ComplementarySubsystemSpaceDimension = ComplementarySubsystemSpace->GetHilbertSpaceDimension();
    }
  cout << "subsystem Hilbert space dimension = " << SubsystemSpaceDimension << endl;

  if (SubsystemSpace != 0)
    {
      for (int i = 0; i < SubsystemSpaceDimension; ++i)
	{
	  unsigned long TmpState = SubsystemSpace->StateDescription[i];
	  double Tmp = 1.0;
	  int TmpPos = 0;
	  while (TmpPos <= SubsystemSpace->LzMax)
	    {
	      switch (TmpState & 03ul)
		{
		case 0x1ul:
		  Tmp *= weightOrbitalADown[TmpPos];
		  break;
		case 0x2ul:
		  Tmp *= weightOrbitalAUp[TmpPos];
		  break;
		case 0x3ul:
		  Tmp *= weightOrbitalADown[TmpPos] * weightOrbitalAUp[TmpPos];
		  break;
		}
	      TmpState >>= 2;
	      ++TmpPos;
	    }
	  for (int j = 0; j < ComplementarySubsystemSpaceDimension; ++j)          
	    entanglementMatrix(i, j) *= Tmp;      
	}
    }

  if (ComplementarySubsystemSpace != 0)
    {
      for (int MinIndex = 0; MinIndex < ComplementarySubsystemSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = ComplementarySubsystemSpace->StateDescription[MinIndex];
	  double Tmp = 1.0;
	  int TmpPos = 0;
	  while (TmpPos <= ComplementarySubsystemSpace->LzMax)
	    {
	      switch (TmpState & 03ul)
		{
		case 0x1ul:
		  Tmp *= weightOrbitalBDown[TmpPos];
		  break;
		case 0x2ul:
		  Tmp *= weightOrbitalBUp[TmpPos];
		  break;
		case 0x3ul:
		  Tmp *= weightOrbitalBDown[TmpPos] * weightOrbitalBUp[TmpPos];
		  break;
		}
	      TmpState >>= 2;
	      ++TmpPos;
	    }
	  for (int j = 0; j < SubsystemSpaceDimension; ++j)
	    entanglementMatrix(j, MinIndex) *= Tmp; 
	}
    }

  if (ComplementarySubsystemSpace != 0)
    {
      delete ComplementarySubsystemSpace;
    }
  if (SubsystemSpace != 0)
    {
      delete SubsystemSpace;
    }
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsystem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// szParitySector = parity sector for the discrete symmetry Sz<->-Sz. Can be either -1 or +1
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& FermionOnSphereWithSpinAllSzSzSymmetry::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector, int szParity,
																       int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
																       int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, 
																       RealMatrix& entanglementMatrix)
{
  int TotalSzSymmetrySector = 1;
  if (this->SzParitySign < 0.0)
    {
      TotalSzSymmetrySector = -1;
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrFermions, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrParticleSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      return entanglementMatrix;	  
    }
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
  
  FermionOnSphereWithSpinAllSzSzSymmetry* SubsystemSpace = 0;
  int SubsystemSpaceDimension = 1;
  if (nbrParticleSector > 0)
    {      
      SubsystemSpace = new FermionOnSphereWithSpinAllSzSzSymmetry(nbrParticleSector, lzSector, nbrOrbitalA - 1, (szParity == -1));
      SubsystemSpaceDimension = SubsystemSpace->GetHilbertSpaceDimension();
    }
  else
    {
      if (szParity == -1)
	{
	  SubsystemSpaceDimension = 0;
	}
    }
  FermionOnSphereWithSpinAllSzSzSymmetry* ComplementarySubsystemSpace = 0;
  int ComplementarySubsystemSpaceDimension = 1;  
  if (ComplementaryNbrParticles > 0)
    {
      ComplementarySubsystemSpace = new FermionOnSphereWithSpinAllSzSzSymmetry(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ((szParity * TotalSzSymmetrySector) == -1));
      ComplementarySubsystemSpaceDimension = ComplementarySubsystemSpace->GetHilbertSpaceDimension();
    }
  else
    {
      if ((szParity * TotalSzSymmetrySector) == -1)
	{
	  ComplementarySubsystemSpaceDimension = 0;
	}
    }
  cout << "subsystem Hilbert space dimension = " << SubsystemSpaceDimension << endl;

  if (SubsystemSpace != 0)
    {
      for (int i = 0; i < SubsystemSpaceDimension; ++i)
	{
	  unsigned long TmpState = SubsystemSpace->StateDescription[i];
	  double Tmp = 1.0;
	  int TmpPos = 0;
	  while (TmpPos <= SubsystemSpace->LzMax)
	    {
	      switch (TmpState & 03ul)
		{
		case 0x1ul:
		  Tmp *= weightOrbitalADown[TmpPos];
		  break;
		case 0x2ul:
		  Tmp *= weightOrbitalAUp[TmpPos];
		  break;
		case 0x3ul:
		  Tmp *= weightOrbitalADown[TmpPos] * weightOrbitalAUp[TmpPos];
		  break;
		}
	      TmpState >>= 2;
	      ++TmpPos;
	    }
	  for (int j = 0; j < ComplementarySubsystemSpaceDimension; ++j)          
	    entanglementMatrix(i, j) *= Tmp;      
	}
    }

  if (ComplementarySubsystemSpace != 0)
    {
      for (int MinIndex = 0; MinIndex < ComplementarySubsystemSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = ComplementarySubsystemSpace->StateDescription[MinIndex];
	  double Tmp = 1.0;
	  int TmpPos = 0;
	  while (TmpPos <= ComplementarySubsystemSpace->LzMax)
	    {
	      switch (TmpState & 03ul)
		{
		case 0x1ul:
		  Tmp *= weightOrbitalBDown[TmpPos];
		  break;
		case 0x2ul:
		  Tmp *= weightOrbitalBUp[TmpPos];
		  break;
		case 0x3ul:
		  Tmp *= weightOrbitalBDown[TmpPos] * weightOrbitalBUp[TmpPos];
		  break;
		}
	      TmpState >>= 2;
	      ++TmpPos;
	    }
	  for (int j = 0; j < SubsystemSpaceDimension; ++j)
	    entanglementMatrix(j, MinIndex) *= Tmp; 
	}
    }

  if (ComplementarySubsystemSpace != 0)
    {
      delete ComplementarySubsystemSpace;
    }
  if (SubsystemSpace != 0)
    {
      delete SubsystemSpace;
    }
  return entanglementMatrix;
}

