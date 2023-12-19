////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on square lattice with spin            //
//                        in filtered momentum space                          //
// using a simplified version more in line with SU(3), SU(6) or SU(12) codes  //
//                                                                            //
//                        last modification : 17/12/2023                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace.h"
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
#include "GeneralTools/MultiColumnASCIIFile.h"

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


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, char* allowedOrbitalsFileName, int kxMomentum, int kyMomentum, unsigned long memory)
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
  this->OrbitalFilteringMask = 0x0ul;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->ParseOrbitalFile(allowedOrbitalsFileName);
  this->LargeHilbertSpaceDimension = this->EvaluateFilteredHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Temporary Hilbert space dimension: " << this->LargeHilbertSpaceDimension << endl;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateFilteredStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      this->FilterHilbertSpace(allowedOrbitalsFileName);
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

// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
// nbrSpinUp = number of particles with spin up
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, char* allowedOrbitalsFileName, int kxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->NbrFermionsUp = nbrSpinUp;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->TotalSpin = (this->NbrFermionsUp - this->NbrFermionsDown);
  this->OrbitalFilteringMask = 0x0ul;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  if (this->NbrFermions <= (this->NbrSiteX * this->NbrSiteY))
    {
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp);
    }
  else
    {
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimensionHoles((2 * this->NbrSiteX * this->NbrSiteY) - this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, (this->NbrSiteX * (this->NbrSiteX - 1)) * this->NbrSiteY, (this->NbrSiteY * (this->NbrSiteY - 1)) * this->NbrSiteX, (this->NbrSiteX * this->NbrSiteY) - this->NbrFermionsUp);      
    }
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Temporary Hilbert space dimension: " << this->LargeHilbertSpaceDimension << endl;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = 0l;
      if (this->NbrFermions <= (this->NbrSiteX * this->NbrSiteY))
	{
	   TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp, 0l);
	}
      else
	{
	  TmpLargeHilbertSpaceDimension = this->GenerateStatesHoles((2 * this->NbrSiteX * this->NbrSiteY) - this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, (this->NbrSiteX * (this->NbrSiteX - 1)) * this->NbrSiteY, (this->NbrSiteY * (this->NbrSiteY - 1)) * this->NbrSiteX, (this->NbrSiteX * this->NbrSiteY) - this->NbrFermionsUp, 0l);
	  SortArrayDownOrdering<unsigned long>(this->StateDescription, TmpLargeHilbertSpaceDimension);
	}
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      this->FilterHilbertSpace(allowedOrbitalsFileName);
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

FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace(const FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace& fermions)
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
  this->SzFlag = fermions.SzFlag;
  this->OrbitalFilteringMask = fermions.OrbitalFilteringMask;
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

FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::~FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace& FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::operator = (const FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace& fermions)
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
  this->OrbitalFilteringMask = fermions.OrbitalFilteringMask;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace(*this);
}

// parse the ascii file providing the orbitals that are allowed
//
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed

void FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::ParseOrbitalFile(char* allowedOrbitalsFileName)
{
  MultiColumnASCIIFile AllowedOrbitalsFile;
  if (AllowedOrbitalsFile.Parse(allowedOrbitalsFileName) == false)
    {
      AllowedOrbitalsFile.DumpErrors(cout) << endl;
      exit(0);
    }
  if (AllowedOrbitalsFile.GetNbrLines() == 0)
    {
      cout << allowedOrbitalsFileName << " is an empty file" << endl;
      exit(0);
    }
  if (AllowedOrbitalsFile.GetNbrColumns() < 3)
    {
      cout << allowedOrbitalsFileName << " needs at least three columns" << endl;
      exit(0);
    }
  int* TmpKxValues = AllowedOrbitalsFile.GetAsIntegerArray(0);
  int* TmpKyValues = AllowedOrbitalsFile.GetAsIntegerArray(1);
  int* TmpBandValues = AllowedOrbitalsFile.GetAsIntegerArray(2);
  
  this->OrbitalFilteringMask = 0x0ul;
  for (int i = 0 ; i < AllowedOrbitalsFile.GetNbrLines(); ++i)
    {
      this->OrbitalFilteringMask |= 0x1ul << ((((TmpKxValues[i]  * this->NbrSiteY) + TmpKyValues[i]) << 1) + TmpBandValues[i]);
    }
  this->OrbitalFilteringMask = ~this->OrbitalFilteringMask;
}

// filter Hilbert to remove forbidden orbitals
//
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed

void FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::FilterHilbertSpace(char* allowedOrbitalsFileName)
{
  this->ParseOrbitalFile(allowedOrbitalsFileName);
  long TmpHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->StateDescription[i] & this->OrbitalFilteringMask) == 0x0ul)
	{
	  TmpHilbertSpaceDimension++;
	}
    }
  unsigned long* TmpStateDescription = new unsigned long [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->StateDescription[i] & this->OrbitalFilteringMask) == 0x0ul)
	{
	  TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	  TmpHilbertSpaceDimension++;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  delete[] this->StateHighestBit;
  this->StateHighestBit = new int [TmpHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::GenerateFilteredStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
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
      unsigned long Mask;
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      Mask = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 1);
	      if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
		{
		  this->StateDescription[pos] = Mask;
		  ++pos;
		}
	      Mask = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 1);;
	      if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
		{
		  this->StateDescription[pos] = Mask;
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
		  Mask = 0x2ul << (((i * this->NbrSiteY) + j) << 1);
		  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
		    {
		      this->StateDescription[pos] = Mask;
		      ++pos;
		    }
		  Mask = 0x1ul << (((i * this->NbrSiteY) + j) << 1);
		  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
		    {
		      this->StateDescription[pos] = Mask;
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }
  unsigned long Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  long TmpPos;
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateFilteredStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
};



// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace::EvaluateFilteredHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      if ((this->OrbitalFilteringMask & (0x2ul << ((((currentKx * this->NbrSiteY) + j) << 1)))) == 0x0ul)
		{
		  Count++;
		}
	      if ((this->OrbitalFilteringMask & (0x1ul << ((((currentKx * this->NbrSiteY) + j) << 1)))) == 0x0ul)
		{
		  Count++;
		}
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  if ((this->OrbitalFilteringMask & (0x2ul << ((((i * this->NbrSiteY) + j) << 1)))) == 0x0ul)
		    {
		      Count++;
		    }
		  if ((this->OrbitalFilteringMask & (0x1ul << ((((i * this->NbrSiteY) + j) << 1)))) == 0x0ul)
		    {
		      Count++;
		    }
		}
	    }
	}
      return Count;
    }
  if ((this->OrbitalFilteringMask & (0x3ul << ((((currentKx * this->NbrSiteY) + currentKy) << 1)))) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
    }
  if ((this->OrbitalFilteringMask & (0x2ul << ((((currentKx * this->NbrSiteY) + currentKy) << 1)))) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
    }
  if ((this->OrbitalFilteringMask & (0x1ul << ((((currentKx * this->NbrSiteY) + currentKy) << 1)))) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
    }
  Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

