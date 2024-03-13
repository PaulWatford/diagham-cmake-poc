////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(3) spin         //
//                          in a filtered momentum space                      //
//                 with a cap on the number of particles per band             //
//                                                                            //
//                        last modification : 02/01/2024                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace.h"
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
#include "GeneralTools/MultiColumnASCIIFile.h"
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


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
// maxNbrParticlesBand0 = maximum number of particles in band 0
// maxNbrParticlesBand1 = maximum number of particles in band 1
// maxNbrParticlesBand2 = maximum number of particles in band 2
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, char* allowedOrbitalsFileName, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, int kxMomentum, int kyMomentum, unsigned long memory)
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
  this->MaxNbrParticlesBand0 = maxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = maxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = maxNbrParticlesBand2;
  this->TotalNbrOrbitals = 3 * this->NbrSiteX * this->NbrSiteY;
  this->MaxTotalMomentumX = 0;
  this->MaxTotalMomentumY = 0;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->OrbitalFilteringMask = 0x0ul;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->ParseOrbitalFile(allowedOrbitalsFileName);
  this->LargeHilbertSpaceDimension = this->EvaluateFilteredHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2);
  if (this->NbrFermions <= (this->TotalNbrOrbitals / 2))
     {
      this->LargeHilbertSpaceDimension = this->EvaluateFilteredHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2);
     }
   else
     {
       int HoleMomentumX = (this->MaxTotalMomentumX - this->KxMomentum) % this->NbrSiteX;
       int HoleMomentumY = (this->MaxTotalMomentumY - this->KyMomentum) % this->NbrSiteY;
       int TmpKx = this->KxMomentum;
       int TmpKy = this->KyMomentum;
       this->KxMomentum = HoleMomentumX;
       this->KyMomentum = HoleMomentumY;
       this->LargeHilbertSpaceDimension = this->EvaluateFilteredHilbertSpaceDimension(this->TotalNbrOrbitals - this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2);
       this->KxMomentum = TmpKx;
       this->KyMomentum = TmpKy;
     }
  cout << "Temporary Hilbert space dimension: " << this->LargeHilbertSpaceDimension << endl;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = 0l;
      if (this->NbrFermions <= (this->TotalNbrOrbitals / 2))
	{
	  TmpLargeHilbertSpaceDimension = this->GenerateFilteredStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2, 0l);
	}
      else
	{
	  int HoleMomentumX = (this->MaxTotalMomentumX - this->KxMomentum) % this->NbrSiteX;
	  int HoleMomentumY = (this->MaxTotalMomentumY - this->KyMomentum) % this->NbrSiteY;
	  int TmpKx = this->KxMomentum;
	  int TmpKy = this->KyMomentum;
	  this->KxMomentum = HoleMomentumX;
	  this->KyMomentum = HoleMomentumY;
	  TmpLargeHilbertSpaceDimension = this->GenerateFilteredStates(this->TotalNbrOrbitals - this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->MaxNbrParticlesBand0, this->MaxNbrParticlesBand1, this->MaxNbrParticlesBand2, 0l);
	  this->KxMomentum = TmpKx;
	  this->KyMomentum = TmpKy;
	  unsigned long TmpMask = ~this->OrbitalFilteringMask;
	  for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	    {
	      this->StateDescription[i] = (~this->StateDescription[i]) & TmpMask;
	    }
	  SortArrayDownOrdering<unsigned long>(this->StateDescription, TmpLargeHilbertSpaceDimension);
	}
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      this->FilterHilbertSpace(allowedOrbitalsFileName);
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace(const FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->MaxNbrParticlesBand0 = fermions.MaxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = fermions.MaxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = fermions.MaxNbrParticlesBand2;
  this->OrbitalFilteringMask = fermions.OrbitalFilteringMask;
  this->TotalNbrOrbitals = fermions.TotalNbrOrbitals;
  this->MaxTotalMomentumX = fermions.MaxTotalMomentumX;
  this->MaxTotalMomentumY = fermions.MaxTotalMomentumY;
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

FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::~FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace& FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::operator = (const FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace& fermions)
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
  this->MaxNbrParticlesBand0 = fermions.MaxNbrParticlesBand0;
  this->MaxNbrParticlesBand1 = fermions.MaxNbrParticlesBand1;
  this->MaxNbrParticlesBand2 = fermions.MaxNbrParticlesBand2;
  this->OrbitalFilteringMask = fermions.OrbitalFilteringMask;
  this->TotalNbrOrbitals = fermions.TotalNbrOrbitals;
  this->MaxTotalMomentumX = fermions.MaxTotalMomentumX;
  this->MaxTotalMomentumY = fermions.MaxTotalMomentumY;
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace(*this);
}

// parse the ascii file providing the orbitals that are allowed
//
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed

void FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::ParseOrbitalFile(char* allowedOrbitalsFileName)
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
  this->TotalNbrOrbitals = AllowedOrbitalsFile.GetNbrLines();
  int TmpMaxNbrParticlesBand0 = 0;
  int TmpMaxNbrParticlesBand1 = 0;
  int TmpMaxNbrParticlesBand2 = 0;
  this->MaxTotalMomentumX = 0;
  this->MaxTotalMomentumY = 0;
  for (int i = 0 ; i < AllowedOrbitalsFile.GetNbrLines(); ++i)
    {
      this->OrbitalFilteringMask |= 0x1ul << ((((TmpKxValues[i]  * this->NbrSiteY) + TmpKyValues[i]) * 3) + TmpBandValues[i]);
      this->MaxTotalMomentumX += TmpKxValues[i];
      this->MaxTotalMomentumY += TmpKyValues[i];
      switch (TmpBandValues[i])
	{
	case 0:
	  TmpMaxNbrParticlesBand0++;
	  break;
	case 1:
	  TmpMaxNbrParticlesBand1++;
	  break;
	case 2:
	  TmpMaxNbrParticlesBand2++;
	  break;
	}
    }
  this->OrbitalFilteringMask = ~this->OrbitalFilteringMask;
  if (this->MaxNbrParticlesBand0 > TmpMaxNbrParticlesBand0)
    {
      this->MaxNbrParticlesBand0 = TmpMaxNbrParticlesBand0;
    }
  if (this->MaxNbrParticlesBand1 > TmpMaxNbrParticlesBand1)
    {
      this->MaxNbrParticlesBand1 = TmpMaxNbrParticlesBand1;
    }
  if (this->MaxNbrParticlesBand2 > TmpMaxNbrParticlesBand2)
    {
      this->MaxNbrParticlesBand2 = TmpMaxNbrParticlesBand2;
    }
}

// filter Hilbert to remove forbidden orbitals
//
// allowedOrbitalsFileName = ascii file providing the orbitals that are allowed

void FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::FilterHilbertSpace(char* allowedOrbitalsFileName)
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
// maxNbrParticlesBand0 = current maximum number of particles in band 0
// maxNbrParticlesBand1 = current maximum number of particles in band 1
// maxNbrParticlesBand2 = current maximum number of particles in band 2
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::GenerateFilteredStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2, long pos)
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
      unsigned long Mask;;
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      Mask = 0x4ul << (((currentKx * this->NbrSiteY) + j) * 3);
	      if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand2 > 0))
		{
		  this->StateDescription[pos] = Mask;
		  ++pos;
		}
	      Mask = 0x2ul << (((currentKx * this->NbrSiteY) + j) * 3);
	      if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand1 > 0))
		{
		  this->StateDescription[pos] = Mask;
		  ++pos;
		}
	      Mask = 0x1ul << (((currentKx * this->NbrSiteY) + j) * 3);
	      if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand0 > 0))
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
		  Mask = 0x4ul << (((i * this->NbrSiteY) + j) * 3);
		  if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand2 > 0))
		    {
		      this->StateDescription[pos] = Mask;
		      ++pos;
		    }
		  Mask = 0x2ul << (((i * this->NbrSiteY) + j) * 3);
		  if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand1 > 0))
		    {
		      this->StateDescription[pos] = Mask;
		      ++pos;
		    }
		  Mask = 0x1ul << (((i * this->NbrSiteY) + j) * 3);
		  if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand0 > 0))
		    {
		      this->StateDescription[pos] = Mask;
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }

  long TmpPos;

  unsigned long Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      TmpPos = this->GenerateFilteredStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }

  return this->GenerateFilteredStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2, pos);
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

long FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace::EvaluateFilteredHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int maxNbrParticlesBand0, int maxNbrParticlesBand1, int maxNbrParticlesBand2)
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
      unsigned long Mask;
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      Mask = 0x4ul << (((currentKx * this->NbrSiteY) + j) * 3);
	      if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand2 > 0))
		{
		  Count++;
		}
	      Mask = 0x2ul << (((currentKx * this->NbrSiteY) + j) * 3);
	      if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand1 > 0))
		{
		  Count++;
		}
	      Mask = 0x1ul << (((currentKx * this->NbrSiteY) + j) * 3);
	      if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand0 > 0))
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
		  Mask = 0x4ul << (((i * this->NbrSiteY) + j) * 3);
		  if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand2 > 0))
		    {
		      Count++;
		    }
		  Mask = 0x2ul << (((i * this->NbrSiteY) + j) * 3);
		  if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand1 > 0))
		    {
		      Count++;
		    }
		  Mask = 0x1ul << (((i * this->NbrSiteY) + j) * 3);
		  if (((this->OrbitalFilteringMask & Mask) == 0x0ul) && (maxNbrParticlesBand0 > 0))
		    {
		      Count++;
		    }
		}
	    }
	}
      return Count;
    }
  unsigned long Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1);
    }
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2 - 1);
    }
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1);
    }
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), maxNbrParticlesBand0 - 1, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2);
    }
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2 - 1);
    }
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0, maxNbrParticlesBand1 - 1, maxNbrParticlesBand2);
    }
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 3);
  if ((this->OrbitalFilteringMask & Mask) == 0x0ul)
    {
      Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, maxNbrParticlesBand0 - 1, maxNbrParticlesBand1, maxNbrParticlesBand2);
    }
  Count += this->EvaluateFilteredHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, maxNbrParticlesBand0, maxNbrParticlesBand1, maxNbrParticlesBand2);
  return Count;
}
