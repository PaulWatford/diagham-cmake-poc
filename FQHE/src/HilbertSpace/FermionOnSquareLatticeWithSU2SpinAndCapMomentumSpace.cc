////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(2) spin         //
//       in momentum space with a cap on the number of particles per band     //
//                                                                            //
//                        last modification : 18/12/2024                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace.h"
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
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "GeneralTools/FilenameTools.h"

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

FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace()
{
}
  
// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// maxNbrParticlesBand0 = maximum number of particles in band 0
// maxNbrParticlesBand1 = maximum number of particles in band 1
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int kxMomentum, int kyMomentum, char* outputDirectory, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->MinNbrParticlesBand0 = 0;
  this->MinNbrParticlesBand1 = 0;
  this->MaxNbrParticlesBand0 = maxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = maxNbrParticlesBand1;
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
  this->TargetSpace = this;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  if (outputDirectory == 0)
    {
      this->GenerateStatesFromSingleBandHilbertSpaces();
    }
  else
    {
      char* TmpName = this->GetDefaultHilbertSpaceFileName();
      char* FullName = new char[strlen(TmpName) + strlen(outputDirectory) + 16];
      sprintf (FullName, "%s/%s", outputDirectory, TmpName);
      if (IsFile(FullName))
	{
	  this->ReadCoreHilbertSpace(FullName);	  
	}
      else
	{
	  this->GenerateStatesFromSingleBandHilbertSpaces();
	  this->WriteCoreHilbertSpace(FullName);
	}
      delete[] TmpName;
      delete[] FullName;
    }
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
// maxNbrParticlesBand0 = maximum number of particles in band 0
// maxNbrParticlesBand1 = maximum number of particles in band 1
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int minNbrParticlesBand0, int minNbrParticlesBand1, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int kxMomentum, int kyMomentum, char* outputDirectory, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->MinNbrParticlesBand0 = minNbrParticlesBand0;
  this->MinNbrParticlesBand1 = minNbrParticlesBand1;
  this->MaxNbrParticlesBand0 = maxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = maxNbrParticlesBand1;
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
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TargetSpace = this;
  if (outputDirectory == 0)
    {
      this->GenerateStatesFromSingleBandHilbertSpaces();
    }
  else
    {
      char* TmpName = this->GetDefaultHilbertSpaceFileName();
      char* FullName = new char[strlen(TmpName) + strlen(outputDirectory) + 16];
      sprintf (FullName, "%s/%s", outputDirectory, TmpName);
      if (IsFile(FullName))
	{
	  this->ReadCoreHilbertSpace(FullName);	  
	}
      else
	{
	  this->GenerateStatesFromSingleBandHilbertSpaces();
	  this->WriteCoreHilbertSpace(FullName);
	}
      delete[] TmpName;
      delete[] FullName;
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


// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace(const FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace& fermions)
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
  this->MaxNbrParticlesBand0 = fermions.MaxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = fermions.MaxNbrParticlesBand1;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
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
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::~FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace& FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::operator = (const FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
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
  this->MaxNbrParticlesBand0 = fermions.MaxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = fermions.MaxNbrParticlesBand1;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace(*this);
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

long FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy)
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

long FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, unsigned long* singleBandStateDescription, long pos)
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
  long TmpPos = this->GenerateSingleBandStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
  unsigned long Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 2);
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

void FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, unsigned long****& singleBandStates)
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
      TmpHoleMask |= 0x1ul << (2 * i);
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

void FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::GenerateStatesFromSingleBandHilbertSpaces()
{
  int TmpMaxBandOccupation = this->MaxNbrParticlesBand0;
  if (TmpMaxBandOccupation < this->MaxNbrParticlesBand1)
    {
      TmpMaxBandOccupation = this->MaxNbrParticlesBand1;
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
      int TmpN1 = this->NbrFermions - TmpN0;
      if ((TmpN1 >= this->MinNbrParticlesBand1) && (TmpN1 <=  this->MaxNbrParticlesBand1))
	{
	  for (int TmpKx0 = 0; TmpKx0 <= TmpSingleBandTotalKxMax[TmpN0]; ++TmpKx0)
	    {
	      for (int TmpKx1 = 0; TmpKx1 <= TmpSingleBandTotalKxMax[TmpN1]; ++TmpKx1)
		{
		  if (((TmpKx0 + TmpKx1) % this->NbrSiteX) == this->KxMomentum)
		    {
		      for (int TmpKy0 = 0; TmpKy0 <= TmpSingleBandTotalKyMax[TmpN0]; ++TmpKy0)
			{
			  if (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0] > 0l)
			    {
			      for (int TmpKy1 = 0; TmpKy1 <= TmpSingleBandTotalKyMax[TmpN1]; ++TmpKy1)
				{
				  if ((((TmpKy0 + TmpKy1) % this->NbrSiteY) == this->KyMomentum) && (TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1] > 0l))
				    {
				      TmpLargeHilbertSpaceDimension += (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0]
									* TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1]);
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
	  int TmpN1 = this->NbrFermions - TmpN0;
	  if ((TmpN1 >= this->MinNbrParticlesBand1) && (TmpN1 <=  this->MaxNbrParticlesBand1))
	    {
	      for (int TmpKx0 = 0; TmpKx0 <= TmpSingleBandTotalKxMax[TmpN0]; ++TmpKx0)
		{
		  for (int TmpKx1 = 0; TmpKx1 <= TmpSingleBandTotalKxMax[TmpN1]; ++TmpKx1)
		    {
		      if (((TmpKx0 + TmpKx1) % this->NbrSiteX) == this->KxMomentum)
			{
			  for (int TmpKy0 = 0; TmpKy0 <= TmpSingleBandTotalKyMax[TmpN0]; ++TmpKy0)
			    {
			      if (TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0] > 0l)
				{
				  for (int TmpKy1 = 0; TmpKy1 <= TmpSingleBandTotalKyMax[TmpN1]; ++TmpKy1)
				    {
				      if ((((TmpKy0 + TmpKy1) % this->NbrSiteY) == this->KyMomentum) && (TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1] > 0l))
					{
					  for (int Pos0 = 0; Pos0 < TmpSingleBandHilbertDimensions[TmpN0][TmpKx0][TmpKy0]; ++Pos0)
					    {
					      for (int Pos1 = 0; Pos1 < TmpSingleBandHilbertDimensions[TmpN1][TmpKx1][TmpKy1]; ++Pos1)
						{
						  this->StateDescription[TmpLargeHilbertSpaceDimension] =  (TmpSingleBandStates[TmpN0][TmpKx0][TmpKy0][Pos0]
													    | (TmpSingleBandStates[TmpN1][TmpKx1][TmpKy1][Pos1] << 1));
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

// provide the default name of the file for Hilbert space storage 
//
// return value = pointer to file name (0 if no default file name exists) 

char* FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace::GetDefaultHilbertSpaceFileName()
{
  char* TmpName = this->FermionOnSquareLatticeWithSU2SpinMomentumSpace::GetDefaultHilbertSpaceFileName();
  if (TmpName == 0)
    {
      return 0;
    }
  char* TmpExtension = new char[256];
  sprintf(TmpExtension, "_minband0_%d_minband1_%d_maxband0_%d_maxband1_%d.hil", this->MinNbrParticlesBand0, this->MinNbrParticlesBand1, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1);
  char* TmpName2 = ReplaceExtensionToFileName(TmpName, ".hil", TmpExtension);
  delete[] TmpExtension;
  delete[] TmpName;
  return TmpName2;
}

