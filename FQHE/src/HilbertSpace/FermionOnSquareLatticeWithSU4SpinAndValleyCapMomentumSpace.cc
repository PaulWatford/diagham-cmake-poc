////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//           class of fermions on a square lattice with SU(4) spin            //
//   in momentum space with a cap on the number of particles in each valley   //
//                                                                            //
//                        last modification : 17/11/2023                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
//

FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace ()
{
  this->MaxNbrParticlesPlus = 0;
  this->MaxNbrParticlesMinus = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// maxNbrParticlesPlus = maximum number of particles in valley plus
// maxNbrParticlesMinus = maximum number of particles in valley minus
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesPlus, int maxNbrParticlesMinus, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->PzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUpPlus = 0;
  this->NbrFermionsDownPlus = 0;
  this->NbrFermionsUpMinus = 0;
  this->NbrFermionsDownMinus = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->MaxNbrParticlesPlus = maxNbrParticlesPlus;
  this->MaxNbrParticlesMinus = maxNbrParticlesMinus;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->HighestBit = (4 * this->NbrSiteX * this->NbrSiteY) - 1;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesPlus, this->MaxNbrParticlesMinus);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l, this->MaxNbrParticlesPlus, this->MaxNbrParticlesMinus);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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
}

// constructor when preserving only spin
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// maxNbrParticlesPlus = maximum number of particles in valley plus
// maxNbrParticlesMinus = maximum number of particles in valley minus
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesPlus, int maxNbrParticlesMinus, int kxMomentum, int kyMomentum, int totalSpin, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->PzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = 0;
  this->NbrFermionsUpPlus = ((this->NbrFermions+this->TotalSpin)/2 + this->TotalIsospin)/2;
  this->NbrFermionsDownPlus = ((this->NbrFermions-this->TotalSpin)/2 + this->TotalIsospin)/2;
  this->NbrFermionsUpMinus = ((this->NbrFermions+this->TotalSpin)/2 - this->TotalIsospin)/2;
  this->NbrFermionsDownMinus = ((this->NbrFermions-this->TotalSpin)/2 - this->TotalIsospin)/2;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->MaxNbrParticlesPlus = maxNbrParticlesPlus;
  this->MaxNbrParticlesMinus = maxNbrParticlesMinus;
  this->HighestBit = (4 * this->NbrSiteX * this->NbrSiteY) - 1;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 (this->NbrFermions + this->TotalSpin) / 2, this->MaxNbrParticlesPlus, this->MaxNbrParticlesMinus);
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l, (this->NbrFermions+this->TotalSpin)/2, this->MaxNbrParticlesPlus, this->MaxNbrParticlesMinus);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << endl;
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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
}


// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace(const FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->SzFlag = fermions.SzFlag;
  this->PzFlag = fermions.PzFlag;
  this->NbrFermionsUpPlus = fermions.NbrFermionsUpPlus;
  this->NbrFermionsDownPlus = fermions.NbrFermionsDownPlus;
  this->NbrFermionsUpMinus = fermions.NbrFermionsUpMinus;
  this->NbrFermionsDownMinus = fermions.NbrFermionsDownMinus;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->HighestBit = fermions.HighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->MaxNbrParticlesPlus = fermions.MaxNbrParticlesPlus;
  this->MaxNbrParticlesMinus = fermions.MaxNbrParticlesMinus;
}

// destructor
//

FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::~FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace& FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::operator = (const FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->PzFlag = fermions.PzFlag;
  this->TotalSpin = fermions.TotalSpin;
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->NbrFermionsUpPlus = fermions.NbrFermionsUpPlus;
  this->NbrFermionsDownPlus = fermions.NbrFermionsDownPlus;
  this->NbrFermionsUpMinus = fermions.NbrFermionsUpMinus;
  this->NbrFermionsDownMinus = fermions.NbrFermionsDownMinus;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->HighestBit = fermions.HighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->MaxNbrParticlesPlus = fermions.MaxNbrParticlesPlus;
  this->MaxNbrParticlesMinus = fermions.MaxNbrParticlesMinus;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace(*this);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// maxNbrParticlesPlus = current maximum number of particles in valley plus
// maxNbrParticlesMinus = current maximum number of particles in valley minus
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int maxNbrParticlesPlus, int maxNbrParticlesMinus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (maxNbrParticlesPlus < 0) || (maxNbrParticlesMinus < 0))
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      if (maxNbrParticlesPlus > 0)
		{// up plus
		  this->StateDescription[pos] = 0x8ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;
		}// up minus
	      if (maxNbrParticlesMinus > 0)
		{
		  this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;
		}
	      if (maxNbrParticlesPlus > 0)
		{// down plus
		  this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;
		}
	      if (maxNbrParticlesMinus > 0)
		{// down minus
		  this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;
		}
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  if (maxNbrParticlesPlus > 0)
		    {
		      this->StateDescription[pos] = 0x8ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		  if (maxNbrParticlesMinus > 0)
		    {
		      this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		  if (maxNbrParticlesPlus > 0)
		    {
		      this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		  if (maxNbrParticlesMinus > 0)
		    {
		      this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }

  // u-p u-m d-p d-m
  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 2);
  unsigned long Mask = 0xful << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
      
  // u-p u-m d-p
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1);
  Mask = 0xeul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p u-m d-m
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2);
  Mask = 0xdul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p u-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0xcul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p d-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1);
  Mask = 0xbul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p d-p
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, maxNbrParticlesPlus, maxNbrParticlesMinus - 2);
  Mask = 0xaul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m d-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2);
  Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m d-p
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m d-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, maxNbrParticlesPlus, maxNbrParticlesMinus - 2);
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, maxNbrParticlesPlus, maxNbrParticlesMinus - 1);
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // d-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // d-p
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, maxNbrParticlesPlus - 1, maxNbrParticlesMinus);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // d-m
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, maxNbrParticlesPlus, maxNbrParticlesMinus - 1);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos, maxNbrParticlesPlus, maxNbrParticlesMinus);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrFermionsUp = current number of fermions with a spin up
// pos = position in StateDescription array where to store states
// maxNbrParticlesPlus = current maximum number of particles in valley plus
// maxNbrParticlesMinus = current maximum number of particles in valley minus
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int nbrFermionsUp, int maxNbrParticlesPlus, int maxNbrParticlesMinus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrFermionsUp < 0) || (nbrFermionsUp > nbrFermions) || (maxNbrParticlesPlus < 0) || (maxNbrParticlesMinus < 0))
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
        {
            this->StateDescription[pos] = 0x0ul;	  
            return (pos + 1l);
        }
      else
	{
          return pos;
	}
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrFermionsUp == 1)
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		  (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  if (maxNbrParticlesPlus > 0)
		    {
		      this->StateDescription[pos] = 0x8ul << (((currentKx * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		  if (maxNbrParticlesMinus > 0)
		    {
		      this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      if (maxNbrParticlesPlus > 0)
			{
			  this->StateDescription[pos] = 0x8ul << (((i * this->NbrSiteY) + j) << 2);
			  ++pos;
			}
		      if (maxNbrParticlesMinus > 0)
			{
			  this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) << 2);
			  ++pos;
			}
		    }
		}
	    }
	}
      else
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		  (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  if (maxNbrParticlesPlus > 0)
		    {
		      this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		  if (maxNbrParticlesMinus > 0)
		    {
		      this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 2);
		      ++pos;
		    }
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      if (maxNbrParticlesPlus > 0)
			{
			  this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 2);
			  ++pos;
			}
		      if (maxNbrParticlesMinus > 0)
			{
			  this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 2);
			  ++pos;		  
			}
		    }
		}
	    }
	}
      return pos;
    }

  // u-p u-m d-p d-m
  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos, nbrFermionsUp - 2, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 2);
  unsigned long Mask = 0xful << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p u-m d-p
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1);
  Mask = 0xeul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p u-m d-m
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2);
  Mask = 0xdul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p u-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 2, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0xcul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p d-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1);
  Mask = 0xbul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p d-p
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, maxNbrParticlesPlus - 2, maxNbrParticlesMinus);
  Mask = 0xaul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  
  // u-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-p
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m d-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2);
  Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m d-p
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  
  // u-m d-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, maxNbrParticlesPlus, maxNbrParticlesMinus - 2);
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // u-m
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1, maxNbrParticlesPlus, maxNbrParticlesMinus - 1);
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // d-p d-m
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  // d-p
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp, maxNbrParticlesPlus - 1, maxNbrParticlesMinus);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  
  // d-m
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp, maxNbrParticlesPlus, maxNbrParticlesMinus - 1);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos, nbrFermionsUp, maxNbrParticlesPlus, maxNbrParticlesMinus);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// maxNbrParticlesPlus = current maximum number of particles in valley plus
// maxNbrParticlesMinus = current maximum number of particles in valley minus
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesPlus, int maxNbrParticlesMinus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (maxNbrParticlesPlus < 0) || (maxNbrParticlesMinus < 0))
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      long TmpIncrement = 0l;
      if (maxNbrParticlesPlus > 0)
	{
	  TmpIncrement += 2l;
	}
      if (maxNbrParticlesMinus > 0)
	{
	  TmpIncrement += 2l;
	}
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += TmpIncrement;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += TmpIncrement;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 2);
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1));
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2));
  
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesPlus - 2, maxNbrParticlesMinus));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesPlus, maxNbrParticlesMinus - 2));
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesPlus - 1, maxNbrParticlesMinus));
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesPlus, maxNbrParticlesMinus - 1));
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, maxNbrParticlesPlus, maxNbrParticlesMinus);
  return Count;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrFermionsUp = current number of fermions with a spin up
// maxNbrParticlesPlus = current maximum number of particles in valley plus
// maxNbrParticlesMinus = current maximum number of particles in valley minus
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrFermionsUp, int maxNbrParticlesPlus, int maxNbrParticlesMinus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)|| (nbrFermionsUp < 0) ||  (nbrFermionsUp > nbrFermions) || (maxNbrParticlesPlus < 0) || (maxNbrParticlesMinus < 0))
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      long TmpIncrement = 0l;
      if (maxNbrParticlesPlus > 0)
	{
	  TmpIncrement++;
	}
      if (maxNbrParticlesMinus > 0)
	{
	  TmpIncrement++;
	}
      
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += TmpIncrement;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += TmpIncrement;
	    }
	}
      return Count;
    }
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrFermionsUp - 2, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 2);
  
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1, maxNbrParticlesPlus - 2, maxNbrParticlesMinus - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 2));
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 2, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, maxNbrParticlesPlus - 2, maxNbrParticlesMinus));
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, maxNbrParticlesPlus, maxNbrParticlesMinus - 2));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp, maxNbrParticlesPlus - 1, maxNbrParticlesMinus - 1);
  
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp, maxNbrParticlesPlus - 1, maxNbrParticlesMinus));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp, maxNbrParticlesPlus, maxNbrParticlesMinus - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1, maxNbrParticlesPlus - 1, maxNbrParticlesMinus));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1, maxNbrParticlesPlus, maxNbrParticlesMinus - 1));
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrFermionsUp, maxNbrParticlesPlus, maxNbrParticlesMinus);
  return Count;
}


