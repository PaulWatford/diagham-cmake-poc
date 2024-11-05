////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(6) spin         //
//     in momentum space  with a cap on the number of particles per band      //
//                                                                            //
//                        last modification : 23/07/2024                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace.h"
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

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace()
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

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSz = 0;
  this->TotalLz = 0;
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
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->GenerateStatesFromSingleBandHilbertSpaces();
  if (this->LargeHilbertSpaceDimension > 0l)
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

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSz = 0;
  this->TotalLz = 0;
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
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->GenerateStatesFromSingleBandHilbertSpaces();
  if (this->LargeHilbertSpaceDimension > 0l)
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

// constructor
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
// totalSz = twice the total Sz (or any U(1) quantum number)
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int minNbrParticlesBand2, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, int totalSz, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSz = totalSz;
  this->TotalLz = 0;
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
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->GenerateSpinConvervedStatesFromSingleBandHilbertSpaces();
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

// copy constructor (without duplicating data)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace(const FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace& fermions)
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
  this->TotalSz = fermions.TotalSz;
  this->SzFlag = fermions.SzFlag;
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

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::~FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace ()
{
}

// assignement (without duplicating data)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace& FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::operator = (const FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace& fermions)
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
  this->TotalSz = fermions.TotalSz;
  this->SzFlag = fermions.SzFlag;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating data)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace(*this);
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

long FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy)
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
  if ((currentKx < 0) || (nbrFermions < 0) || (((nbrFermions * currentKx) + currentTotalKx) < singleBandTotalKx)
      || (((nbrFermions * (this->NbrSiteY - 1)) + currentTotalKy) < singleBandTotalKy))
    return 0l;
  long Count = 0;
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), singleBandTotalKx, singleBandTotalKy);
  Count += (2 * this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, singleBandTotalKx, singleBandTotalKy));
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, singleBandTotalKx, singleBandTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for a single band
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalSz = current total spin 
// singleBandTotalKx = total momentum along x
// singleBandTotalKy = total momentum along y
// singleBandTotalSz = total spin 
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentTotalSz, int singleBandTotalKx, int singleBandTotalKy, int singleBandTotalSz)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == singleBandTotalKx) && (currentTotalKy == singleBandTotalKy) && (currentTotalSz == singleBandTotalSz))
	return 1l;
      else	
	return 0l;
    }
  if ((currentKx < 0) || (nbrFermions < 0) || (((nbrFermions * currentKx) + currentTotalKx) < singleBandTotalKx)
      || (((nbrFermions * (this->NbrSiteY - 1)) + currentTotalKy) < singleBandTotalKy))
    return 0l;
  long Count = 0;
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalSz, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz);
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalSz + 1, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz);
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalSz - 1, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz);  
  Count += this->EvaluateSingleBandHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentTotalSz, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz);
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

long FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, unsigned long* singleBandStateDescription, long pos)
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
  if ((currentKx < 0) || (nbrFermions < 0) || (((nbrFermions * currentKx) + currentTotalKx) < singleBandTotalKx)
      || (((nbrFermions * (this->NbrSiteY - 1)) + currentTotalKy) < singleBandTotalKy))
    return pos;
  long TmpPos = this->GenerateSingleBandStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
  unsigned long Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) * 6);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;

  TmpPos = this->GenerateSingleBandStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) * 6);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;

  TmpPos = this->GenerateSingleBandStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 6);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;

  return this->GenerateSingleBandStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
}

// generate all states corresponding to the constraints for a single band
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalSz = current total spin 
// singleBandTotalKx = total momentum along x
// singleBandTotalKy = total momentum along y
// singleBandTotalSz = total spin 
// singleBandStateDescription = pointer to the single band state description array
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentTotalSz, int singleBandTotalKx, int singleBandTotalKy, int singleBandTotalSz, unsigned long* singleBandStateDescription, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == singleBandTotalKx) && (currentTotalKy == singleBandTotalKy) && (currentTotalSz == singleBandTotalSz))
	{
	  singleBandStateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if ((currentKx < 0) || (nbrFermions < 0) || (((nbrFermions * currentKx) + currentTotalKx) < singleBandTotalKx)
      || (((nbrFermions * (this->NbrSiteY - 1)) + currentTotalKy) < singleBandTotalKy))
    return pos;
  long TmpPos = this->GenerateSingleBandStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalSz, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz, singleBandStateDescription, pos);
  unsigned long Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) * 6);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;

  TmpPos = this->GenerateSingleBandStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalSz + 1, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz, singleBandStateDescription, pos);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) * 6);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;

  TmpPos = this->GenerateSingleBandStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalSz - 1, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz, singleBandStateDescription, pos);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 6);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;

  return this->GenerateSingleBandStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentTotalSz, singleBandTotalKx, singleBandTotalKy, singleBandTotalSz, singleBandStateDescription, pos);
}

// evaluate all the single band Hilbert spaces
//
// maxBandOccupation = maiximum occupation of a single band
// singleBandTotalKxMax = reference on the array for the maximum total Kx values 
// singleBandTotalKyMax = reference on the array for the maximum total Ky values
// singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky)
// singleBandStates = reference on the array for all the Hilbert space basis states 

void FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, unsigned long****& singleBandStates)
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
      TmpHoleMask |= 0x9ul << (3 * i);
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

// evaluate all the single band Hilbert spaces
//
// maxBandOccupation = maximum occupation of a single band
// singleBandTotalKxMax = reference on the array for the maximum total Kx values 
// singleBandTotalKyMax = reference on the array for the maximum total Ky values
// singleBandTotalSz = reference on the array for the maximum total Sz values
// singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky, fourth index being the total Sz)
// singleBandStates = reference on the array for all the Hilbert space basis states 

void FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, int*& singleBandTotalSz, long****& singleBandHilbertDimensions, unsigned long*****& singleBandStates)
{
  singleBandTotalKxMax = new int[maxBandOccupation + 1];
  singleBandTotalKyMax = new int[maxBandOccupation + 1];
  singleBandTotalSz = new int[maxBandOccupation + 1];
  singleBandHilbertDimensions = new long***[maxBandOccupation + 1];
  singleBandStates = new unsigned long****[maxBandOccupation + 1];
  int TmpHalfMaxBandOccupation = (this->NbrSiteX * this->NbrSiteY);
  if (maxBandOccupation < TmpHalfMaxBandOccupation)
    {
      TmpHalfMaxBandOccupation = maxBandOccupation;
    }
  // needed to deactivate construction for ph conjugation 
  //  TmpHalfMaxBandOccupation = maxBandOccupation;
  for (int i = 0; i <= TmpHalfMaxBandOccupation; ++i)
    {
      singleBandTotalKxMax[i] = (this->NbrSiteX - 1) * i;
      singleBandTotalKyMax[i] = (this->NbrSiteY - 1) * i;
      singleBandTotalSz[i] = i;
      singleBandHilbertDimensions[i] = new long**[singleBandTotalKxMax[i] + 1];
      singleBandStates[i] = new unsigned long***[singleBandTotalKxMax[i] + 1];
      for (int j = 0; j <= singleBandTotalKxMax[i]; ++j)
	{
	  singleBandHilbertDimensions[i][j] = new long*[singleBandTotalKyMax[i] + 1];
	  singleBandStates[i][j] = new unsigned long**[singleBandTotalKyMax[i] + 1];
	  for (int k = 0; k <= singleBandTotalKyMax[i]; ++k)
	    {
	      singleBandHilbertDimensions[i][j][k] = new long[singleBandTotalSz[i] + 1];
	      singleBandStates[i][j][k] = new unsigned long*[singleBandTotalSz[i] + 1];
	      for (int l = 0; l <= singleBandTotalSz[i]; ++l)
		{
		  singleBandHilbertDimensions[i][j][k][l] = this->EvaluateSingleBandHilbertSpaceDimension(i, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0, j, k, (2 * l) - i);
		  //		  cout << "N=" << i << " Kx=" << j << " Ky=" << k << " 2Sz=" << ((2 * l) - i) << " dim=" << singleBandHilbertDimensions[i][j][k][l] << " | " <<  singleBandTotalKxMax[i] << " " << singleBandTotalKyMax[i] << " " << singleBandTotalSz[i] << endl;
		  if (singleBandHilbertDimensions[i][j][k][l] > 0l)
		    {
		      singleBandStates[i][j][k][l] = new unsigned long[singleBandHilbertDimensions[i][j][k][l]];
		      long Tmp = this->GenerateSingleBandStates(i, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0, j, k, (2 * l) - i, singleBandStates[i][j][k][l], 0l);
		      if (Tmp != singleBandHilbertDimensions[i][j][k][l])
			{
			  cout << "Single band Hilbert space generation error at N=" << i << " Kx=" << j << " Ky=" << k << " 2Sz=" << ((2 * l) - i) << endl;
			}
		    }
		  else
		    {
		      singleBandStates[i][j][k][l] = 0;
		    }
		}
	    }
	}
    }
  
  unsigned long TmpHoleMask = 0x0ul;
  for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
    {
      TmpHoleMask |= 0x9ul << (6 * i);
    }

  for (int i = TmpHalfMaxBandOccupation + 1; i <= maxBandOccupation; ++i)
    {
      singleBandTotalKxMax[i] = (this->NbrSiteX - 1) * i;
      singleBandTotalKyMax[i] = (this->NbrSiteY - 1) * i;
      singleBandTotalSz[i] = i;
      singleBandHilbertDimensions[i] = new long**[singleBandTotalKxMax[i] + 1];
      singleBandStates[i] = new unsigned long***[singleBandTotalKxMax[i] + 1];
      for (int j = 0; j <= singleBandTotalKxMax[i]; ++j)
	{
	  singleBandHilbertDimensions[i][j] = new long*[singleBandTotalKyMax[i] + 1];
	  singleBandStates[i][j] = new unsigned long**[singleBandTotalKyMax[i] + 1];
	  for (int k = 0; k <= singleBandTotalKyMax[i]; ++k)
	    {
	      singleBandHilbertDimensions[i][j][k] = new long[singleBandTotalSz[i] + 1];
	      singleBandStates[i][j][k] = new unsigned long*[singleBandTotalSz[i] + 1];
	      for (int l = 0; l <= singleBandTotalSz[i]; ++l)
		{
		  int TmpNbrHoles = (2 * this->NbrSiteX * this->NbrSiteY) - i;
		  int TmpHolesTotalKx = (2 * this->NbrSiteY * ((this->NbrSiteX * (this->NbrSiteX - 1)) >> 1)) - j;
		  int TmpHolesTotalKy = (2 * this->NbrSiteX * ((this->NbrSiteY * (this->NbrSiteY - 1)) >> 1)) - k;
		  int TmpHolesTotalSz = (TmpNbrHoles - ((2 * l) - i)) / 2;
		  // cout << "Nh=" << TmpNbrHoles << " kxh=" << TmpHolesTotalKx << " kyh=" << TmpHolesTotalKy << " szh=" << TmpHolesTotalSz << " | " << "N=" << i << " kx=" << j << " ky=" << k << " sz=" << l << " | " << singleBandTotalKxMax[TmpNbrHoles] << " " << singleBandTotalKyMax[TmpNbrHoles] << " " << singleBandTotalSz[TmpNbrHoles] << endl;
		   if ((singleBandTotalKxMax[TmpNbrHoles] < TmpHolesTotalKx) || (singleBandTotalKyMax[TmpNbrHoles] < TmpHolesTotalKy)
		      || (TmpHolesTotalKx < 0) || (TmpHolesTotalKy < 0)
		      || (TmpNbrHoles < TmpHolesTotalSz) || (TmpHolesTotalSz < 0) || (singleBandTotalSz[TmpNbrHoles] < TmpHolesTotalSz))
		    {
		      singleBandHilbertDimensions[i][j][k][l] = 0l;
		      singleBandStates[i][j][k][l] = 0;
		    }
		  else
		    {
		      singleBandHilbertDimensions[i][j][k][l] = singleBandHilbertDimensions[TmpNbrHoles][TmpHolesTotalKx][TmpHolesTotalKy][TmpHolesTotalSz];
		      //		      cout << "test " << singleBandHilbertDimensions[i][j][k][l] << endl;
		      if (singleBandHilbertDimensions[i][j][k][l] > 0l)
			{
			  long Tmp = singleBandHilbertDimensions[i][j][k][l];
			  singleBandStates[i][j][k][l] = new unsigned long[Tmp];
			  unsigned long* TmpInput = singleBandStates[TmpNbrHoles][TmpHolesTotalKx][TmpHolesTotalKy][TmpHolesTotalSz];
			  unsigned long* TmpOutput = singleBandStates[i][j][k][l];
			  for (int m = 0; m < Tmp; ++m)
			    {
			      TmpOutput[m] = (~TmpInput[m]) & TmpHoleMask;
			    }
			}
		      else
			{
			  singleBandStates[i][j][k][l] = 0;
			}
		    }
		}
	    }
	}
    }
}

// generate all states using the single band hilbert spaces
//

void FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::GenerateStatesFromSingleBandHilbertSpaces()
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
  if (TmpMaxBandOccupation > (2 * this->NbrSiteX * this->NbrSiteY))
    {
      TmpMaxBandOccupation = 2 * this->NbrSiteX * this->NbrSiteY;
    }
   if (TmpMaxBandOccupation > NbrFermions)
    {
      TmpMaxBandOccupation = NbrFermions;
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


// generate all states using the single band hilbert spaces for the spin conserved case
//

void FermionOnSquareLatticeWithSU6SpinAndCapMomentumSpace::GenerateSpinConvervedStatesFromSingleBandHilbertSpaces()
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
  TmpMaxBandOccupation *= 2; 
  if (TmpMaxBandOccupation > (2 * this->NbrSiteX * this->NbrSiteY))
    {
      TmpMaxBandOccupation = 2 * this->NbrSiteX * this->NbrSiteY;
    }
  if (TmpMaxBandOccupation > NbrFermions)
    {
      TmpMaxBandOccupation = NbrFermions;
    }

  int TmpNbrSpinUp = (this->NbrFermions + this->TotalSz) / 2;
  
  int* TmpSingleBandTotalKxMax = 0;
  int* TmpSingleBandTotalKyMax = 0;
  int* TmpSingleBandTotalSz = 0;
  long**** TmpSingleBandHilbertDimensions = 0;
  unsigned long***** TmpSingleBandStates = 0;
  this->GenerateAllSingleBandHilbertSpaces(TmpMaxBandOccupation, TmpSingleBandTotalKxMax, TmpSingleBandTotalKyMax, TmpSingleBandTotalSz, TmpSingleBandHilbertDimensions, TmpSingleBandStates);

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
				  for (int TmpKy1 = 0; TmpKy1 <= TmpSingleBandTotalKyMax[TmpN1]; ++TmpKy1)
				    {
				      for (int TmpKy2 = 0; TmpKy2 <= TmpSingleBandTotalKyMax[TmpN2]; ++TmpKy2)
					{
					  if (((TmpKy0 + TmpKy1 + TmpKy2) % this->NbrSiteY) == this->KyMomentum)
					    {
					      for (int TmpSz0 = 0; TmpSz0 <= TmpSingleBandTotalSz[TmpN0]; ++TmpSz0)
						{
						  if (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0][TmpSz0] > 0l)
						    {
						      for (int TmpSz1 = 0; TmpSz1 <= TmpSingleBandTotalSz[TmpN1]; ++TmpSz1)
							{
							  if (TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1][TmpSz1] > 0l)
							    {
							      for (int TmpSz2 = 0; TmpSz2 <= TmpSingleBandTotalSz[TmpN2]; ++TmpSz2)
								{								  
								  if (((TmpSz0 + TmpSz1 + TmpSz2) == TmpNbrSpinUp) && (TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2][TmpSz2] > 0l))
								    {
								      TmpLargeHilbertSpaceDimension += (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0][TmpSz0]
													* TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1][TmpSz1]
													* TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2][TmpSz2]);
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
				      for (int TmpKy1 = 0; TmpKy1 <= TmpSingleBandTotalKyMax[TmpN1]; ++TmpKy1)
					{
					  for (int TmpKy2 = 0; TmpKy2 <= TmpSingleBandTotalKyMax[TmpN2]; ++TmpKy2)
					    {
					      if (((TmpKy0 + TmpKy1 + TmpKy2) % this->NbrSiteY) == this->KyMomentum)
						{
						  for (int TmpSz0 = 0; TmpSz0 <= TmpSingleBandTotalSz[TmpN0]; ++TmpSz0)
						    {
						      if (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0][TmpSz0] > 0l)
							{
							  for (int TmpSz1 = 0; TmpSz1 <= TmpSingleBandTotalSz[TmpN1]; ++TmpSz1)
							    {
							      if (TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1][TmpSz1] > 0l)
								{
								  for (int TmpSz2 = 0; TmpSz2 <= TmpSingleBandTotalSz[TmpN2]; ++TmpSz2)
								    {
								      if (((TmpSz0 + TmpSz1 + TmpSz2) == TmpNbrSpinUp) && (TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2][TmpSz2] > 0l))
									{
									  for (int Pos0 = 0; Pos0 < TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0][TmpSz0]; ++Pos0)
									    {
									      for (int Pos1 = 0; Pos1 < TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1][TmpSz1]; ++Pos1)
										{
										  for (int Pos2 = 0; Pos2 < TmpSingleBandHilbertDimensions[TmpN2][TmpKx2][TmpKy2][TmpSz2]; ++Pos2)
										    {
										      this->StateDescription[TmpLargeHilbertSpaceDimension] =  (TmpSingleBandStates[TmpN0][TmpKx0][TmpKy0][TmpSz0][Pos0]
																		| (TmpSingleBandStates[TmpN1][TmpKx1][TmpKy1][TmpSz1][Pos1] << 1)
																		| (TmpSingleBandStates[TmpN2][TmpKx2][TmpKy2][TmpSz2][Pos2] << 2));
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
	      for (int l = 0; l <= TmpSingleBandTotalSz[i]; ++l)
		{
		  if (TmpSingleBandStates[i][j][k][l] != 0)
		    {
		      delete[] TmpSingleBandStates[i][j][k][l];
		    }
		}
	      delete[] TmpSingleBandHilbertDimensions[i][j][k];
	      delete[] TmpSingleBandStates[i][j][k];
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
