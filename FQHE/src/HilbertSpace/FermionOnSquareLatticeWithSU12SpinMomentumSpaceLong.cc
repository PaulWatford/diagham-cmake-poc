////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(12) spin        //
//                 in momentum space and for more than 5 orbitals            //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong.h"
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


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSz = 0;
  this->PzFlag = false;
  this->TotalPz = 0;
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// totalSz = twice the total Sz (or any U(1) quantum number)
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSz, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSz = totalSz;
  this->PzFlag = false;
  this->TotalPz = 0;
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 (this->NbrFermions + this->TotalSz) / 2);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
								(this->NbrFermions + this->TotalSz) / 2, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// totalSz = twice the total Sz (or any U(1) quantum number)
// totalPz = twice the total Pz (or any U(1) quantum number)
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSz, int totalPz, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSz = totalSz;
  this->PzFlag = true;
  this->TotalPz = totalPz;
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 (this->NbrFermions + this->TotalPz) / 2, (this->NbrFermions + this->TotalSz) / 2);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
								(this->NbrFermions + this->TotalPz) / 2, (this->NbrFermions + this->TotalSz) / 2, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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

FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong(const FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong& fermions)
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
  this->TotalSz = fermions.TotalSz;
  this->SzFlag = fermions.SzFlag;
  this->PzFlag = fermions.PzFlag;
  this->TotalPz = fermions.TotalPz;
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

FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::~FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong ()
{
}

// assignement (without duplicating data)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong& FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::operator = (const FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong& fermions)
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
  this->TotalSz = fermions.TotalSz;
  this->PzFlag = fermions.PzFlag;
  this->TotalPz = fermions.TotalPz;
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::Clone()
{
  return new FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i * 12));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & ((ULONGLONG) 0x800ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Au+)";
      if ((Tmp & ((ULONGLONG) 0x400ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Bu+)";
      if ((Tmp & ((ULONGLONG) 0x200ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Cu+)";
      if ((Tmp & ((ULONGLONG) 0x100ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Au-)";
      if ((Tmp & ((ULONGLONG) 0x80ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Bu-)";
      if ((Tmp & ((ULONGLONG) 0x40ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Cu-)";
      if ((Tmp & ((ULONGLONG) 0x20ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Ad+)";
      if ((Tmp & ((ULONGLONG) 0x10ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Bd+)";
      if ((Tmp & ((ULONGLONG) 0x8ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Cd+)";
      if ((Tmp & ((ULONGLONG) 0x4ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Ad-)";
      if ((Tmp & ((ULONGLONG) 0x2ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Bd-)";
      if ((Tmp & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",Cd-)";
    }
  Str << "]";
  return Str;
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

long FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
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
	      int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
	      this->StateDescription[pos] = ((ULONGLONG) 0x800ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x400ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x200ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x100ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x80ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x40ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x20ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x10ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << TmpShift;
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << TmpShift;
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = ((ULONGLONG) 0x800ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x400ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x200ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x100ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x80ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x40ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x20ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x10ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << TmpShift;
		  ++pos;
		}
	    }
	}
      return pos;
    }


  long TmpPos = 0l;

  for (int TmpMask = 0xfff; TmpMask > 0; --TmpMask)
    {
      int TmpNbrParticles = TmpMask & 1;
      TmpNbrParticles += (TmpMask >> 1) & 1;
      TmpNbrParticles += (TmpMask >> 2) & 1;
      TmpNbrParticles += (TmpMask >> 3) & 1;
      TmpNbrParticles += (TmpMask >> 4) & 1;
      TmpNbrParticles += (TmpMask >> 5) & 1;
      TmpNbrParticles += (TmpMask >> 6) & 1;
      TmpNbrParticles += (TmpMask >> 7) & 1;
      TmpNbrParticles += (TmpMask >> 8) & 1;
      TmpNbrParticles += (TmpMask >> 9) & 1;
      TmpNbrParticles += (TmpMask >> 10) & 1;
      TmpNbrParticles += (TmpMask >> 11) & 1;
      if (TmpNbrParticles <= nbrFermions)
	{
	  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
	  ULONGLONG Mask = ((ULONGLONG) TmpMask) << (((currentKx * this->NbrSiteY) + currentKy) * 12);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
    }

  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
};


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of particles with a spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrSpinUp < 0) || (nbrSpinUp > nbrFermions))
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrSpinUp == 1)
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = ((ULONGLONG) 0x800ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x400ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x200ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x100ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x80ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x40ul) << TmpShift;
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		      this->StateDescription[pos] = ((ULONGLONG) 0x800ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x400ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x200ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x100ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x80ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x40ul) << TmpShift;
		      ++pos;
		    }
		}
	    }
	}
      else
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = 0x800ul << TmpShift;
		  this->StateDescription[pos] = ((ULONGLONG) 0x20ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x10ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << TmpShift;
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		      this->StateDescription[pos] = ((ULONGLONG) 0x20ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x10ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << TmpShift;
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }


  long TmpPos = 0l;

  for (int TmpMask = 0xfff; TmpMask > 0; --TmpMask)
    {
      int TmpNbrParticles = (TmpMask >> 11) & 1;
      TmpNbrParticles += (TmpMask >> 10) & 1;
      TmpNbrParticles += (TmpMask >> 9) & 1;
      TmpNbrParticles += (TmpMask >> 8) & 1;
      TmpNbrParticles += (TmpMask >> 7) & 1;
      TmpNbrParticles += (TmpMask >> 6) & 1;
      int TmpNbrParticlesUp = TmpNbrParticles;
      TmpNbrParticles += (TmpMask >> 5) & 1;
      TmpNbrParticles += (TmpMask >> 4) & 1;
      TmpNbrParticles += (TmpMask >> 3) & 1;
      TmpNbrParticles += (TmpMask >> 2) & 1;
      TmpNbrParticles += (TmpMask >> 1) & 1;
      TmpNbrParticles += TmpMask & 1;
      if ((TmpNbrParticles <= nbrFermions) && (TmpNbrParticlesUp <= nbrSpinUp))
	{
	  TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1, currentTotalKx + (TmpNbrParticles * currentKx), currentTotalKy + (TmpNbrParticles * currentKy), nbrSpinUp - TmpNbrParticlesUp, pos);
	  ULONGLONG Mask = ((ULONGLONG) TmpMask) << (((currentKx * this->NbrSiteY) + currentKy) * 12);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescription[pos] |= Mask;
	    }
	}
    }

  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrSpinUp, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrValleyPlus = number of particles in valley plus
// nbrSpinUp = number of particles with a spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrValleyPlus, int nbrSpinUp, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrSpinUp < 0) || (nbrSpinUp > nbrFermions) || (nbrValleyPlus < 0) || (nbrValleyPlus > nbrFermions))
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
      if ((nbrSpinUp == 1) && (nbrValleyPlus == 1))
	{
	   for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = ((ULONGLONG) 0x800ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x400ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x200ul) << TmpShift;
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		      this->StateDescription[pos] = ((ULONGLONG) 0x800ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x400ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x200ul) << TmpShift;
		      ++pos;
		    }
		}
	    }
	  return pos;
	}
      if ((nbrSpinUp == 1) && (nbrValleyPlus == 0))
	{
	   for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = ((ULONGLONG) 0x100ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x80ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x40ul) << TmpShift;
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		      this->StateDescription[pos] = ((ULONGLONG) 0x100ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x80ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x40ul) << TmpShift;
		      ++pos;
		    }
		}
	    }
	  return pos;
	}
       if ((nbrSpinUp == 0) && (nbrValleyPlus == 1))
	 {     
	   for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = ((ULONGLONG) 0x20ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x10ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << TmpShift;
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		      this->StateDescription[pos] = ((ULONGLONG) 0x20ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x10ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << TmpShift;
		      ++pos;
		    }
		}
	    }
	  return pos;
	}
      else
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  int TmpShift = (((currentKx * this->NbrSiteY) + j) * 12);
		  this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << TmpShift;
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << TmpShift;
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      int TmpShift = (((i * this->NbrSiteY) + j) * 12);
		      this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << TmpShift;
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << TmpShift;
		      ++pos;
		    }
		}
	    }
	  return pos;
	}
    }


  long TmpPos = 0l;

  for (int TmpMask = 0xfff; TmpMask > 0; --TmpMask)
    {
      int TmpNbrParticles = (TmpMask >> 11) & 1;
      TmpNbrParticles += (TmpMask >> 10) & 1;
      TmpNbrParticles += (TmpMask >> 9) & 1;
      int TmpNbrParticlesPlus = TmpNbrParticles;
      TmpNbrParticles += (TmpMask >> 8) & 1;
      TmpNbrParticles += (TmpMask >> 7) & 1;
      TmpNbrParticles += (TmpMask >> 6) & 1;
      int TmpNbrParticlesUp = TmpNbrParticles;
      int TmpNbrParticles2 = (TmpMask >> 5) & 1;
      TmpNbrParticles2 += (TmpMask >> 4) & 1;
      TmpNbrParticles2 += (TmpMask >> 3) & 1;
      TmpNbrParticlesPlus += TmpNbrParticles2;
      TmpNbrParticles += TmpNbrParticles2;
      TmpNbrParticles += (TmpMask >> 2) & 1;
      TmpNbrParticles += (TmpMask >> 1) & 1;
      TmpNbrParticles += TmpMask & 1;
      if ((TmpNbrParticles <= nbrFermions) && (TmpNbrParticlesUp <= nbrSpinUp) && (TmpNbrParticlesPlus <= nbrValleyPlus))
	{
	  TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1, currentTotalKx + (TmpNbrParticles * currentKx), currentTotalKy + (TmpNbrParticles * currentKy), nbrValleyPlus - TmpNbrParticlesPlus, nbrSpinUp - TmpNbrParticlesUp, pos);
	  ULONGLONG Mask = ((ULONGLONG) TmpMask) << (((currentKx * this->NbrSiteY) + currentKy) * 12);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescription[pos] |= Mask;
	    }
	}
    }

  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrValleyPlus, nbrSpinUp, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
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
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += 12l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 12l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 12, currentKx, currentKy - 1, currentTotalKx + (12 * currentKx), currentTotalKy + (12 * currentKy));
  Count += (12 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy)));
  Count += (66 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy)));
  Count += (220 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy)));
  Count += (495 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy)));
  Count += (792 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy)));
  Count += (924 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy)));
  Count += (792 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy)));
  Count += (495 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy)));
  Count += (220 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (66 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (12 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy)));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of particles with a spin up
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrSpinUp < 0) || (nbrSpinUp > nbrFermions))
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
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += 6l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 6l;
	    }
	}
      return Count;
    }
  if (nbrFermions >= 12)
    {
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 12, currentKx, currentKy - 1, currentTotalKx + (12 * currentKx), currentTotalKy + (12 * currentKy), nbrSpinUp - 6));
    }

  if (nbrFermions >= 11)
    {
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy), nbrSpinUp - 6));
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy), nbrSpinUp - 5));
    }

  if (nbrFermions >= 10)
    {
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrSpinUp - 6));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrSpinUp - 5));
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrSpinUp - 4));
    }

  if (nbrFermions >= 9)
    {
      Count += (20 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrSpinUp - 6));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrSpinUp - 5));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrSpinUp - 4));
      Count += (20 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrSpinUp - 3));
    }

  if (nbrFermions >= 8)
    {
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrSpinUp - 6));
      Count += (120 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrSpinUp - 5));
      Count += (225 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrSpinUp - 4));
      Count += (120 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrSpinUp - 3));
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrSpinUp - 2));
    }

  if (nbrFermions >= 7)
    {
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrSpinUp - 6));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrSpinUp - 5));
      Count += (300 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrSpinUp - 4));
      Count += (300 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrSpinUp - 3));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrSpinUp - 2));
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrSpinUp - 1));
    }

  if (nbrFermions >= 6)
    {
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp - 6));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp - 5));
      Count += (225 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp - 4));
      Count += (400 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp - 3));
      Count += (225 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp - 1));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrSpinUp));
    }

  if (nbrFermions >= 5)
    {
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrSpinUp - 5));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrSpinUp - 4));
      Count += (300 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrSpinUp - 3));
      Count += (300 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrSpinUp - 2));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrSpinUp - 1));
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrSpinUp));
    }

  if (nbrFermions >= 4)
    {
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrSpinUp - 4));
      Count += (120 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrSpinUp - 3));
      Count += (225 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrSpinUp - 2));
      Count += (120 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrSpinUp - 1));
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrSpinUp));
    }

  if (nbrFermions >= 3)
    {
      Count += (20 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrSpinUp - 3));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrSpinUp - 2));
      Count += (90 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrSpinUp - 1));
      Count += (20 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrSpinUp));
    }

  if (nbrFermions >= 2)
    {
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrSpinUp - 1));
      Count += (15 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrSpinUp));
    }

  if (nbrFermions >= 1)
    {
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy), nbrSpinUp - 1));
      Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy), nbrSpinUp));
    }

  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrSpinUp);
  return Count;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrValleyPlus = number of particles in valley plus 
// nbrSpinUp = number of particles with a spin up
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrValleyPlus, int nbrSpinUp)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrSpinUp < 0) || (nbrSpinUp > nbrFermions) || (nbrValleyPlus < 0) || (nbrValleyPlus > nbrFermions))
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
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += 3l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 3l;
	    }
	}
      return Count;
    }
  
  if (nbrFermions >= 12)
    {
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 12, currentKx, currentKy - 1, currentTotalKx + (12 * currentKx), currentTotalKy + (12 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 6));
    }

  if (nbrFermions >= 11)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 6));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 6));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 5));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 11, currentKx, currentKy - 1, currentTotalKx + (11 * currentKx), currentTotalKy + (11 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 5));
    }

  if (nbrFermions >= 10)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 6));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 5));
      Count += (18 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 5));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 4));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 10, currentKx, currentKy - 1, currentTotalKx + (10 * currentKx), currentTotalKy + (10 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 4));
    }

  if (nbrFermions >= 9)
    {
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 6));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 5));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 5));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 4));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 4));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 4));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 4));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 3));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 9, currentKx, currentKy - 1, currentTotalKx + (9 * currentKx), currentTotalKy + (9 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
    }

  if (nbrFermions >= 8)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 6));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 6));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 5));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 5));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 5));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 5));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 4));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 4));
      Count += (99 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 4));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 4));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 3));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 3));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 3));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 3));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 2));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 2));
    }

  if (nbrFermions >= 7)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 6));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 5));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 5));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 5));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 4));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 4));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 4));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 4));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 3));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 3));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 3));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 3));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 2));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 1));
    }

  if (nbrFermions >= 6)
    {
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 6));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 5));
      Count += (18 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 4));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 4));
      Count += (99 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 4));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 4));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 4));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 6, nbrSpinUp - 3));
      Count += (18 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 3));
      Count += (99 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 3));
      Count += (164 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
      Count += (99 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 3));
      Count += (18 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 3));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 2));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 2));
      Count += (99 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 2));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 1));
      Count += (18 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 1));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 1));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy), nbrValleyPlus - 3, nbrSpinUp));
    }

  if (nbrFermions >= 5)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 5));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 5));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 4));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 4));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 4));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 3));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 3));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 3));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 3));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus, nbrSpinUp - 3));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 5, nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 2));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 2));
      Count += (111 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 2));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 1));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 1));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 1));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 3, nbrSpinUp));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy), nbrValleyPlus - 2, nbrSpinUp));
    }

  if (nbrFermions >= 4)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 4));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 4));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 3));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 3));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 3));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 2));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 2));
      Count += (99 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 2));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus, nbrSpinUp - 2));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 4, nbrSpinUp - 1));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 1));
      Count += (54 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 1));
      Count += (30 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 3, nbrSpinUp));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 2, nbrSpinUp));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrValleyPlus - 1, nbrSpinUp));
    }

  if (nbrFermions >= 3)
    {
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 3));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus, nbrSpinUp - 3));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 2));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 3, nbrSpinUp - 1));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 1));
      Count += (36 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 1));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus, nbrSpinUp - 1));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 3, nbrSpinUp));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 2, nbrSpinUp));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus - 1, nbrSpinUp));
      Count += (1 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrValleyPlus, nbrSpinUp));
    }

  if (nbrFermions >= 2)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 2));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus, nbrSpinUp - 2));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus - 2, nbrSpinUp - 1));
      Count += (18 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 1));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus - 2, nbrSpinUp));
      Count += (9 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus - 1, nbrSpinUp));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrValleyPlus, nbrSpinUp));
    }

  if (nbrFermions >= 1)
    {
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy), nbrValleyPlus - 1, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy), nbrValleyPlus, nbrSpinUp - 1));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy), nbrValleyPlus - 1, nbrSpinUp));
      Count += (3 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + (1 * currentKx), currentTotalKy + (1 * currentKy), nbrValleyPlus, nbrSpinUp));
    }

  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrValleyPlus, nbrSpinUp);
  return Count;
}

