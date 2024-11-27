////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(3) spin         //
//       in momentum space with a cap on the number of particles per band     //
//                                                                            //
//                        last modification : 03/12/2023                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace.h"
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
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

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

FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace()
{
  this->MinNbrParticlesBand0 = 0;
  this->MinNbrParticlesBand1 = 0;
  this->MinNbrParticlesBand2 = 0;

  this->MaxNbrParticlesBand0 = 0;
  this->MaxNbrParticlesBand1 = 0;
  this->MaxNbrParticlesBand2 = 0;
}
  
// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// maxNbrParticlesBand0 = maximum number of particles in band 0
// maxNbrParticlesBand1 = maximum number of particles in band 1
// maxNbrParticlesBand2 = maximum number of particles in band 2
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->PzFlag = false;
  this->TotalLz = 0;
  this->TotalTz = 0;
  this->TotalY = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->MinNbrParticlesBand0 = 0;
  this->MinNbrParticlesBand1 = 0;
  this->MinNbrParticlesBand2 = 0;
  this->MaxNbrParticlesBand0 = maxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = maxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = maxNbrParticlesBand2;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  if (this->MaxNbrParticlesBand0 > this->LzMax)
    {
      this->MaxNbrParticlesBand0 = this->LzMax;
    }
  if (this->MaxNbrParticlesBand1 > this->LzMax)
    {
      this->MaxNbrParticlesBand1 = this->LzMax;
    }
  if (this->MaxNbrParticlesBand2 > this->LzMax)
    {
      this->MaxNbrParticlesBand2 = this->LzMax;
    }
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  if (true)
    {
      this->GenerateStatesFromSingleBandHilbertSpaces();
    }
  else
    {
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2);
      if (this->LargeHilbertSpaceDimension >= (1l << 30))
	this->HilbertSpaceDimension = 0;
      else
	this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
	  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
	  long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2, 0l);
	  if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	    {
	      cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	    }
	}
    }
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
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

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// minNbrParticlesBand0 = minimum number of particles in band 0
// minNbrParticlesBand1 = minimum number of particles in band 1
// minNbrParticlesBand2 = minimum number of particles in band 2
// maxNbrParticlesBand0 = maximum number of particles in band 0
// maxNbrParticlesBand1 = maximum number of particles in band 1
// maxNbrParticlesBand2 = maximum number of particles in band 2
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->PzFlag = false;
  this->TotalLz = 0;
  this->TotalTz = 0;
  this->TotalY = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->MinNbrParticlesBand0 = minNbrParticlesBand0;
  this->MinNbrParticlesBand1 = minNbrParticlesBand1;
  this->MinNbrParticlesBand2 = minNbrParticlesBand2;
  this->MaxNbrParticlesBand0 = maxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = maxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = maxNbrParticlesBand2;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  if (this->MaxNbrParticlesBand0 > this->LzMax)
    {
      this->MaxNbrParticlesBand0 = this->LzMax;
    }
  if (this->MaxNbrParticlesBand1 > this->LzMax)
    {
      this->MaxNbrParticlesBand1 = this->LzMax;
    }
  if (this->MaxNbrParticlesBand2 > this->LzMax)
    {
      this->MaxNbrParticlesBand2 = this->LzMax;
    }
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->GenerateStatesFromSingleBandHilbertSpaces();
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
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

FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace(const FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->MinNbrParticlesBand0 = fermions.MinNbrParticlesBand0;
  this->MinNbrParticlesBand1 = fermions.MinNbrParticlesBand1;
  this->MinNbrParticlesBand2 = fermions.MinNbrParticlesBand2;
  this->MaxNbrParticlesBand0 = fermions.MaxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = fermions.MaxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = fermions.MaxNbrParticlesBand2;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::~FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace& FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::operator = (const FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace& fermions)
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
  this->MinNbrParticlesBand0 = fermions.MinNbrParticlesBand0;
  this->MinNbrParticlesBand1 = fermions.MinNbrParticlesBand1;
  this->MinNbrParticlesBand2 = fermions.MinNbrParticlesBand2;
  this->MaxNbrParticlesBand0 = fermions.MaxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = fermions.MaxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = fermions.MaxNbrParticlesBand2;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->PzFlag = fermions.PzFlag;
  this->TotalTz = fermions.TotalTz;
  this->TotalY  = fermions.TotalY;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// maxNbrParticlesBand0 = current maximum number of particles in band 0
// maxNbrParticlesBand1 = current maximum number of particles in band 1
// maxNbrParticlesBand2 = current maximum number of particles in band 2
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (maxNbrParticlesBand0 < 0) || (maxNbrParticlesBand1 < 0) || (maxNbrParticlesBand2 < 0))
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
	      if (maxNbrParticlesBand2 > 0)
		{
		  this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) * 3);
		  ++pos;
		}
	      if (maxNbrParticlesBand1 > 0)
		{
		  this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) * 3);
		  ++pos;
		}
	      if (maxNbrParticlesBand0 > 0)
		{
		  this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) * 3);
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
		  if (maxNbrParticlesBand2 > 0)
		    {
		      this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) * 3);
		      ++pos;
		    }
		  if (maxNbrParticlesBand1 > 0)
		    {
		      this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) * 3);
		      ++pos;
		    }
		  if (maxNbrParticlesBand0 > 0)
		    {
		      this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) * 3);
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }


  long TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1, pos);
  unsigned long Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1, pos);
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1, pos);
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1, pos);
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2, pos);
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2, pos);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2, pos);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// maxNbrParticlesBand0 = current maximum number of particles in band 0
// maxNbrParticlesBand1 = current maximum number of particles in band 1
// maxNbrParticlesBand2 = current maximum number of particles in band 2
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (maxNbrParticlesBand0 < 0) || (maxNbrParticlesBand1 < 0) || (maxNbrParticlesBand2 < 0))
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
      if (maxNbrParticlesBand0 > 0)
	{
	  TmpIncrement++;
	}
      if (maxNbrParticlesBand1 > 0)
	{
	  TmpIncrement++;
	}
      if (maxNbrParticlesBand2 > 0)
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
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1);

  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2));

  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2));
  Count += (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2));

  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2);
  return Count;
}

// evaluate Hilbert space dimension for a single band
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// singleBandTotalKx = total momentum along x
// singleBandTotalKy = total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == singleBandTotalKx) && (currentTotalKy == singleBandTotalKy))
	return 1l;
      else	
	return 0l;
    }
  if ((currentKx < 0) || (((nbrFermions * currentKx) + currentTotalKx) < singleBandTotalKx)
      || (((nbrFermions * (this->NbrSiteY - 1)) + currentTotalKy) < singleBandTotalKy))
    return 0l;
  long Count = 0;
  // if (nbrFermions == 1)
  //   {
  //     if (((currentKx + currentTotalKx) == singleBandTotalKx) && ((currentKy + currentTotalKy) == singleBandTotalKy))
  // 	++Count;
  //     return Count;
  //   }
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, singleBandTotalKx, singleBandTotalKy);
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, singleBandTotalKx, singleBandTotalKy);
  return Count;
}


// generate all states corresponding to the constraints for a single band
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// singleBandTotalKx = total momentum along x
// singleBandTotalKy = total momentum along y
// singleBandStateDescription = pointer to the single band state description array
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, unsigned long* singleBandStateDescription, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == singleBandTotalKx) && (currentTotalKy == singleBandTotalKy))
	{
	  singleBandStateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  // if (nbrFermions == 1)
  //   {
  //     if (((currentKx + currentTotalKx) == singleBandTotalKx) && ((currentKy + currentTotalKy) == singleBandTotalKy))
  // 	{
  // 	  singleBandStateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  // 	  ++pos;
  // 	}
  //     return pos;
  //   }
  long TmpPos = this->GenerateSingleBandStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
  unsigned long Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;
  return this->GenerateSingleBandStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
};

// evaluate all the single band Hilbert spaces
//
// maxBandOccupation = maiximum occupation of a single band
// singleBandTotalKxMax = reference on the array for the maximum total Kx values 
// singleBandTotalKyMax = reference on the array for the maximum total Ky values
// singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky)
// singleBandStates = reference on the array for all the Hilbert space basis states 

void FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, unsigned long****& singleBandStates)
{
  singleBandTotalKxMax = new int[maxBandOccupation + 1];
  singleBandTotalKyMax = new int[maxBandOccupation + 1];
  singleBandHilbertDimensions = new long**[maxBandOccupation + 1];
  singleBandStates = new unsigned long***[maxBandOccupation + 1];
  int TmpHalfMaxBandOccupation = (this->NbrSiteX * this->NbrSiteY) >> 1;
  if (maxBandOccupation < TmpHalfMaxBandOccupation)
    {
      TmpHalfMaxBandOccupation = maxBandOccupation;
    }
  for (int i = 0; i <= TmpHalfMaxBandOccupation; ++i)
    //  for (int i = 0; i <= maxBandOccupation; ++i)
    {
      singleBandTotalKxMax[i] = (this->NbrSiteX - 1) * i;
      singleBandTotalKyMax[i] = (this->NbrSiteY - 1) * i;
      singleBandHilbertDimensions[i] = new long*[singleBandTotalKxMax[i] + 1];
      singleBandStates[i] = new unsigned long**[singleBandTotalKxMax[i] + 1];
      for (int j = 0; j <= singleBandTotalKxMax[i]; ++j)
	{
	  singleBandHilbertDimensions[i][j] = new long[singleBandTotalKyMax[i] + 1];
	  singleBandStates[i][j] = new unsigned long*[singleBandTotalKyMax[i] + 1];
	  for (int k = 0; k <= singleBandTotalKyMax[i]; ++k)
	    {
	      singleBandHilbertDimensions[i][j][k] = this->EvaluateSingleBandHilbertSpaceDimension(i, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, j, k);
	      if (singleBandHilbertDimensions[i][j][k] > 0l)
		{
		  singleBandStates[i][j][k] = new unsigned long[singleBandHilbertDimensions[i][j][k]];
		  long Tmp = this->GenerateSingleBandStates(i, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, j, k, singleBandStates[i][j][k], 0l);
		  if (Tmp != singleBandHilbertDimensions[i][j][k])
		    {
		      cout << "Single band Hilbert space generation error at N=" << i << " Kx=" << j << " Ky=" << k << endl;
		    }
		}
	      else
		{
		  singleBandStates[i][j][k] = 0;
		}
	    }
	}
    }
  
  unsigned long TmpHoleMask = 0x0ul;
  for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
    {
      TmpHoleMask |= 0x1ul << (3 * i);
    }
  
  for (int i = TmpHalfMaxBandOccupation + 1; i <= maxBandOccupation; ++i)
    {
      singleBandTotalKxMax[i] = (this->NbrSiteX - 1) * i;
      singleBandTotalKyMax[i] = (this->NbrSiteY - 1) * i;
      singleBandHilbertDimensions[i] = new long*[singleBandTotalKxMax[i] + 1];
      singleBandStates[i] = new unsigned long**[singleBandTotalKxMax[i] + 1];
      for (int j = 0; j <= singleBandTotalKxMax[i]; ++j)
	{
	  singleBandHilbertDimensions[i][j] = new long[singleBandTotalKyMax[i] + 1];
	  singleBandStates[i][j] = new unsigned long*[singleBandTotalKyMax[i] + 1];
	  for (int k = 0; k <= singleBandTotalKyMax[i]; ++k)
	    {
	      int TmpNbrHoles = (this->NbrSiteX * this->NbrSiteY) - i;
	      int TmpHolesTotalKx = (this->NbrSiteY * ((this->NbrSiteX * (this->NbrSiteX - 1)) >> 1)) - j;
	      int TmpHolesTotalKy = (this->NbrSiteX * ((this->NbrSiteY * (this->NbrSiteY - 1)) >> 1)) - k;
	      if ((singleBandTotalKxMax[TmpNbrHoles] < TmpHolesTotalKx) || (singleBandTotalKyMax[TmpNbrHoles] < TmpHolesTotalKy)
		  || (TmpHolesTotalKx < 0) || (TmpHolesTotalKy < 0))
		{
		  singleBandHilbertDimensions[i][j][k] = 0l;
		  singleBandStates[i][j][k] = 0;
		}
	      else
		{
		  singleBandHilbertDimensions[i][j][k] = singleBandHilbertDimensions[TmpNbrHoles][TmpHolesTotalKx][TmpHolesTotalKy];
		  if (singleBandHilbertDimensions[i][j][k] > 0l)
		    {
		      long Tmp = singleBandHilbertDimensions[i][j][k];
		      singleBandStates[i][j][k] = new unsigned long[Tmp];
		      unsigned long* TmpInput = singleBandStates[TmpNbrHoles][TmpHolesTotalKx][TmpHolesTotalKy];
		      unsigned long* TmpOutput = singleBandStates[i][j][k];
		      for (int l = 0; l < Tmp; ++l)
			{
			  TmpOutput[l] = (~TmpInput[l]) & TmpHoleMask;
			}
		    }
		  else
		    {
		      singleBandStates[i][j][k] = 0;
		    }
		}
	    }
	}
    }
}

// generate all states using the single band hilbert spaces
//

void FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::GenerateStatesFromSingleBandHilbertSpaces()
{
  int TmpMaxBandOccupation = this->MaxNbrParticlesBand0;
  if (TmpMaxBandOccupation < this->MaxNbrParticlesBand1)
    {
      TmpMaxBandOccupation = this->MaxNbrParticlesBand1;
    }
  if (TmpMaxBandOccupation < this->MaxNbrParticlesBand2)
    {
      TmpMaxBandOccupation = this->MaxNbrParticlesBand2;
    }
  if (TmpMaxBandOccupation > NbrFermions)
    {
      TmpMaxBandOccupation = NbrFermions;
    }
  if (TmpMaxBandOccupation > (this->NbrSiteX * this->NbrSiteY))
    {
      TmpMaxBandOccupation = this->NbrSiteX * this->NbrSiteY;
    }
  
  int* TmpSingleBandTotalKxMax = 0;
  int* TmpSingleBandTotalKyMax = 0;
  long*** TmpSingleBandHilbertDimensions = 0;
  unsigned long**** TmpSingleBandStates = 0;
  this->GenerateAllSingleBandHilbertSpaces(TmpMaxBandOccupation, TmpSingleBandTotalKxMax, TmpSingleBandTotalKyMax, TmpSingleBandHilbertDimensions, TmpSingleBandStates);
  
  long TmpLargeHilbertSpaceDimension = 0l;
  for (int TmpN0 = this->MinNbrParticlesBand0; TmpN0 <= this->MaxNbrParticlesBand0; ++TmpN0)
    {
      for (int TmpN1 = this->MinNbrParticlesBand1; TmpN1 <= this->MaxNbrParticlesBand1; ++TmpN1)
	{
	  int TmpN2 = this->NbrFermions - TmpN0 - TmpN1;
	  if ((TmpN2 >= this->MinNbrParticlesBand2) && (TmpN2 <=  this->MaxNbrParticlesBand2))
	    {
	      for (int TmpKx0 = 0; TmpKx0 <= TmpSingleBandTotalKxMax[TmpN0]; ++TmpKx0)
		{
		  for (int TmpKx1 = 0; TmpKx1 <= TmpSingleBandTotalKxMax[TmpN1]; ++TmpKx1)
		    {
		      for (int TmpKx2 = 0; TmpKx2 <= TmpSingleBandTotalKxMax[TmpN2]; ++TmpKx2)
			{
			  if (((TmpKx0 + TmpKx1 + TmpKx2) % this->NbrSiteX) == this->KxMomentum)
			    {
			      for (int TmpKy0 = 0; TmpKy0 <= TmpSingleBandTotalKyMax[TmpN0]; ++TmpKy0)
				{
				  if (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0] > 0l)
				    {
				      for (int TmpKy1 = 0; TmpKy1 <= TmpSingleBandTotalKyMax[TmpN1]; ++TmpKy1)
					{
					  if (TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1] > 0l)
					    {
					      for (int TmpKy2 = 0; TmpKy2 <= TmpSingleBandTotalKyMax[TmpN2]; ++TmpKy2)
						{
						  if ((((TmpKy0 + TmpKy1 + TmpKy2) % this->NbrSiteY) == this->KyMomentum) && (TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2] > 0l))
						    {
						      //						      cout << TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0] << " " <<  TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1] << " " << TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2] << endl;
						      TmpLargeHilbertSpaceDimension += (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0]
											* TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1]
											* TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2]);
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;  
  cout << "Temporary Hilbert space dimension=" << TmpLargeHilbertSpaceDimension << endl;
  if (TmpLargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
      TmpLargeHilbertSpaceDimension = 0l;
      for (int TmpN0 = this->MinNbrParticlesBand0; TmpN0 <= this->MaxNbrParticlesBand0; ++TmpN0)
	{
	  for (int TmpN1 = this->MinNbrParticlesBand1; TmpN1 <= this->MaxNbrParticlesBand1; ++TmpN1)
	    {
	      int TmpN2 = this->NbrFermions - TmpN0 - TmpN1;
	      if ((TmpN2 >= this->MinNbrParticlesBand2) && (TmpN2 <=  this->MaxNbrParticlesBand2))
		{
		  for (int TmpKx0 = 0; TmpKx0 <= TmpSingleBandTotalKxMax[TmpN0]; ++TmpKx0)
		    {
		      for (int TmpKx1 = 0; TmpKx1 <= TmpSingleBandTotalKxMax[TmpN1]; ++TmpKx1)
			{
			  for (int TmpKx2 = 0; TmpKx2 <= TmpSingleBandTotalKxMax[TmpN2]; ++TmpKx2)
			    {
			      if (((TmpKx0 + TmpKx1 + TmpKx2) % this->NbrSiteX) == this->KxMomentum)
				{
				  for (int TmpKy0 = 0; TmpKy0 <= TmpSingleBandTotalKyMax[TmpN0]; ++TmpKy0)
				    {
				      if (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0] > 0l)
					{
					  for (int TmpKy1 = 0; TmpKy1 <= TmpSingleBandTotalKyMax[TmpN1]; ++TmpKy1)
					    {
					      if (TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1] > 0l)
						{
						  for (int TmpKy2 = 0; TmpKy2 <= TmpSingleBandTotalKyMax[TmpN2]; ++TmpKy2)
						    {
						      if ((((TmpKy0 + TmpKy1 + TmpKy2) % this->NbrSiteY) == this->KyMomentum) && (TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2] > 0l))
							{
							  for (int Pos0 = 0; Pos0 < TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0]; ++Pos0)
							    {
							      for (int Pos1 = 0; Pos1 < TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1]; ++Pos1)
								{
								  for (int Pos2 = 0; Pos2 < TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2]; ++Pos2)
								    {
								      this->StateDescription[TmpLargeHilbertSpaceDimension] =  (TmpSingleBandStates[TmpN0][TmpKx0][TmpKy0][Pos0]
																| (TmpSingleBandStates[TmpN1][TmpKx1][TmpKy1][Pos1] << 1)
																| (TmpSingleBandStates[TmpN2][TmpKx2][TmpKy2][Pos2] << 2));
								      TmpLargeHilbertSpaceDimension++;
								    }
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
      SortArrayDownOrdering<unsigned long>(this->StateDescription, TmpLargeHilbertSpaceDimension);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      else
	{
	  cout << "Hilbert space dimension " << this->LargeHilbertSpaceDimension << endl;
	}
      if (this->LargeHilbertSpaceDimension >= (1l << 30))
	this->HilbertSpaceDimension = 0;
      else
	this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
     }
  for (int i = 0; i <= TmpMaxBandOccupation; ++i)
    {
      for (int j = 0; j <= TmpSingleBandTotalKxMax[i]; ++j)
	{
	  for (int k = 0; k <= TmpSingleBandTotalKyMax[i]; ++k)
	    {
	      if (TmpSingleBandStates[i][j][k] != 0)
		{
		  delete[] TmpSingleBandStates[i][j][k];
		}
	    }
	  delete[] TmpSingleBandHilbertDimensions[i][j];
	  delete[] TmpSingleBandStates[i][j];
	}
      delete[] TmpSingleBandHilbertDimensions[i];
      delete[] TmpSingleBandStates[i];
    }
  delete[] TmpSingleBandHilbertDimensions;
  delete[] TmpSingleBandStates;
  delete[] TmpSingleBandTotalKxMax;
  delete[] TmpSingleBandTotalKyMax;
}    


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;

  int SubsystemMaxNbrParticlesBand0 = this->MaxNbrParticlesBand0;
  int SubsystemMaxNbrParticlesBand1 = this->MaxNbrParticlesBand1;
  int SubsystemMaxNbrParticlesBand2 = this->MaxNbrParticlesBand2;
  if (SubsystemMaxNbrParticlesBand0 > nbrParticleSector)
    {
      SubsystemMaxNbrParticlesBand0 = nbrParticleSector;
    }
  if (SubsystemMaxNbrParticlesBand1 > nbrParticleSector)
    {
      SubsystemMaxNbrParticlesBand1 = nbrParticleSector;
    }
  if (SubsystemMaxNbrParticlesBand2 > nbrParticleSector)
    {
      SubsystemMaxNbrParticlesBand2 = nbrParticleSector;
    }
  int ComplementaryMaxNbrParticlesBand0 = this->MaxNbrParticlesBand0;
  int ComplementaryMaxNbrParticlesBand1 = this->MaxNbrParticlesBand1;
  int ComplementaryMaxNbrParticlesBand2 = this->MaxNbrParticlesBand2;
  if (ComplementaryMaxNbrParticlesBand0 > ComplementaryNbrParticles)
    {
      ComplementaryMaxNbrParticlesBand0 = ComplementaryNbrParticles;
    }
  if (ComplementaryMaxNbrParticlesBand1 > ComplementaryNbrParticles)
    {
      ComplementaryMaxNbrParticlesBand1 = ComplementaryNbrParticles;
    }
  if (ComplementaryMaxNbrParticlesBand2 > ComplementaryNbrParticles)
    {
      ComplementaryMaxNbrParticlesBand2 = ComplementaryNbrParticles;
    }

  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  cout << "maxband0 = " << this->MaxNbrParticlesBand0 << " " << SubsystemMaxNbrParticlesBand0 << " " << ComplementaryMaxNbrParticlesBand0 << endl;
  cout << "maxband1 = " << this->MaxNbrParticlesBand1 << " " << SubsystemMaxNbrParticlesBand1 << " " << ComplementaryMaxNbrParticlesBand1 << endl;
  cout << "maxband2 = " << this->MaxNbrParticlesBand2 << " " << SubsystemMaxNbrParticlesBand2 << " " << ComplementaryMaxNbrParticlesBand2 << endl;
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, SubsystemMaxNbrParticlesBand0, SubsystemMaxNbrParticlesBand1, SubsystemMaxNbrParticlesBand2, kxSector, kySector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementaryMaxNbrParticlesBand0, ComplementaryMaxNbrParticlesBand1, ComplementaryMaxNbrParticlesBand2, ComplementaryKxMomentum, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}
  
// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
												    int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 0.0;
	  for (int i = 0; i < nbrGroundStates; ++i)
	    TmpDensityMatrix(0, 0) += weights[i];
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 0.0;
	  for (int i = 0; i < nbrGroundStates; ++i)
	    TmpDensityMatrix(0, 0) += weights[i];
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;

  int SubsystemMaxNbrParticlesBand0 = this->MaxNbrParticlesBand0;
  int SubsystemMaxNbrParticlesBand1 = this->MaxNbrParticlesBand1;
  int SubsystemMaxNbrParticlesBand2 = this->MaxNbrParticlesBand2;
  if (SubsystemMaxNbrParticlesBand0 > nbrParticleSector)
    {
      SubsystemMaxNbrParticlesBand0 = nbrParticleSector;
    }
  if (SubsystemMaxNbrParticlesBand1 > nbrParticleSector)
    {
      SubsystemMaxNbrParticlesBand1 = nbrParticleSector;
    }
  if (SubsystemMaxNbrParticlesBand2 > nbrParticleSector)
    {
      SubsystemMaxNbrParticlesBand2 = nbrParticleSector;
    }
  int ComplementaryMaxNbrParticlesBand0 = this->MaxNbrParticlesBand0;
  int ComplementaryMaxNbrParticlesBand1 = this->MaxNbrParticlesBand1;
  int ComplementaryMaxNbrParticlesBand2 = this->MaxNbrParticlesBand2;
  if (ComplementaryMaxNbrParticlesBand0 > ComplementaryNbrParticles)
    {
      ComplementaryMaxNbrParticlesBand0 = ComplementaryNbrParticles;
    }
  if (ComplementaryMaxNbrParticlesBand1 > ComplementaryNbrParticles)
    {
      ComplementaryMaxNbrParticlesBand1 = ComplementaryNbrParticles;
    }
  if (ComplementaryMaxNbrParticlesBand2 > ComplementaryNbrParticles)
    {
      ComplementaryMaxNbrParticlesBand2 = ComplementaryNbrParticles;
    }
  
  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  cout << "maxband0 = " << this->MaxNbrParticlesBand0 << " " << SubsystemMaxNbrParticlesBand0 << " " << ComplementaryMaxNbrParticlesBand0 << endl;
  cout << "maxband1 = " << this->MaxNbrParticlesBand1 << " " << SubsystemMaxNbrParticlesBand1 << " " << ComplementaryMaxNbrParticlesBand1 << endl;
  cout << "maxband2 = " << this->MaxNbrParticlesBand2 << " " << SubsystemMaxNbrParticlesBand2 << " " << ComplementaryMaxNbrParticlesBand2 << endl;
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, SubsystemMaxNbrParticlesBand0, SubsystemMaxNbrParticlesBand1, SubsystemMaxNbrParticlesBand2, kxSector, kySector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementaryMaxNbrParticlesBand0, ComplementaryMaxNbrParticlesBand1, ComplementaryMaxNbrParticlesBand2, ComplementaryKxMomentum, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, nbrGroundStates, groundStates, weights, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

