////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(4) spin         //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 26/09/2011                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
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

FermionOnSquareLatticeWithSU4SpinMomentumSpace::FermionOnSquareLatticeWithSU4SpinMomentumSpace ()
{
  this->NbrSiteX = 0;
  this->NbrSiteY = 0;
  this->KxMomentum = 0;
  this->KyMomentum = 0;
  this->SzFlag = false;
  this->PzFlag = false;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinMomentumSpace::FermionOnSquareLatticeWithSU4SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
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
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->HighestBit = (4 * this->NbrSiteX * this->NbrSiteY) - 1;
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
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
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
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinMomentumSpace::FermionOnSquareLatticeWithSU4SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, unsigned long memory)
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
  this->HighestBit = (4 * this->NbrSiteX * this->NbrSiteY) - 1;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 (this->NbrFermions + this->TotalSpin) / 2);
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l, (this->NbrFermions+this->TotalSpin)/2);
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


// constructor when conserving spin and isospin
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinMomentumSpace::FermionOnSquareLatticeWithSU4SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->PzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->NbrFermionsUpPlus = ((this->NbrFermions+this->TotalSpin)/2 + this->TotalIsospin)/2;
  this->NbrFermionsDownPlus = ((this->NbrFermions-this->TotalSpin)/2 + this->TotalIsospin)/2;
  this->NbrFermionsUpMinus = ((this->NbrFermions+this->TotalSpin)/2 - this->TotalIsospin)/2;
  this->NbrFermionsDownMinus = ((this->NbrFermions-this->TotalSpin)/2 - this->TotalIsospin)/2;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->HighestBit = (4 * this->NbrSiteX * this->NbrSiteY) - 1;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, (this->NbrFermions+this->TotalSpin)/2, (this->NbrFermions+this->TotalIsospin)/2);
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l, (this->NbrFermions+this->TotalSpin)/2, (this->NbrFermions+this->TotalIsospin)/2);
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

// constructor when preserving the three Cartan quantum numbers
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// totalEntanglement = twice the total entanglement value
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinMomentumSpace::FermionOnSquareLatticeWithSU4SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin,
												int totalEntanglement, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->PzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->TotalEntanglement = totalEntanglement;
  this->NbrFermionsUpPlus = (this->NbrFermions + this->TotalSpin + this->TotalIsospin + this->TotalEntanglement) / 4;
  this->NbrFermionsUpMinus = (this->NbrFermions + this->TotalSpin - this->TotalIsospin - this->TotalEntanglement) / 4;
  this->NbrFermionsDownPlus = (this->NbrFermions - this->TotalSpin + this->TotalIsospin - this->TotalEntanglement) / 4;
  this->NbrFermionsDownMinus = (this->NbrFermions - this->TotalSpin - this->TotalIsospin + this->TotalEntanglement) / 4;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->HighestBit = (4 * this->NbrSiteX * this->NbrSiteY) - 1;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 this->NbrFermionsDownMinus, this->NbrFermionsDownPlus, this->NbrFermionsUpMinus, this->NbrFermionsUpPlus);
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
								this->NbrFermionsDownMinus, this->NbrFermionsDownPlus, this->NbrFermionsUpMinus, this->NbrFermionsUpPlus, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
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

FermionOnSquareLatticeWithSU4SpinMomentumSpace::FermionOnSquareLatticeWithSU4SpinMomentumSpace(const FermionOnSquareLatticeWithSU4SpinMomentumSpace& fermions)
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
}

// destructor
//

FermionOnSquareLatticeWithSU4SpinMomentumSpace::~FermionOnSquareLatticeWithSU4SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU4SpinMomentumSpace& FermionOnSquareLatticeWithSU4SpinMomentumSpace::operator = (const FermionOnSquareLatticeWithSU4SpinMomentumSpace& fermions)
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
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeWithSU4SpinMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU4SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeWithSU4SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 2));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & 0x8ul) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",u,A)";
      if ((Tmp & 0x4ul) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",u,B)";
      if ((Tmp & 0x2ul) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",d,A)";
      if ((Tmp & 0x1ul) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",d,B)";
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

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
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
	      this->StateDescription[pos] = 0x8ul << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	      this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	      this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	      this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = 0x8ul << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		}
	    }
	}
      return pos;
    }

  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos);
  unsigned long Mask = 0xful << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = 0xeul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = 0xdul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = 0xcul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = 0xbul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = 0xaul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
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
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int nbrFermionsUp)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrFermionsUp < 0) || (nbrFermionsUp > nbrFermions))
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
		  this->StateDescription[pos] = 0x8ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;		  
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      this->StateDescription[pos] = 0x8ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;
		      this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;		  
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
		  this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 2);
		  ++pos;		  
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		    {
		      this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;
		      this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 2);
		      ++pos;		  
		    }
		}
	    }
	}
      return pos;
    }

  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos, nbrFermionsUp - 2);
  unsigned long Mask = 0xful << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2);
  Mask = 0xeul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2);
  Mask = 0xdul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 2);
  Mask = 0xcul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1);
  Mask = 0xbul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1);
  Mask = 0xaul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1);
  Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1);
  Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1);
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1);
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1);
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp);
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos, nbrFermionsUp);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrFermionsUp = current number of fermions with a spin up
// nbrFermionsPlus = current number of fermions with a plus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int nbrFermionsUp, int nbrFermionsPlus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrFermionsUp < 0) || (nbrFermionsPlus < 0) ||  
      (nbrFermionsUp > nbrFermions) || (nbrFermionsPlus > nbrFermions))
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
            switch (nbrFermionsUp + 2*nbrFermionsPlus)
            {
                case 3:
                    this->StateDescription[pos] = 0x8ul << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
                case 1:
                    this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
                case 2:
                    this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
                case 0:
                    this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
            }
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
            {
                switch (nbrFermionsUp + 2*nbrFermionsPlus)
                {
                    case 3:
                        this->StateDescription[pos] = 0x8ul << (((i * this->NbrSiteY) + j) << 2);
                        break;
                    case 1:
                        this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) << 2);
                        break;
                    case 2:
                        this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 2);
                        break;
                    case 0:
                        this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 2);
                        break;
                }
                ++pos;
            }
	    }
	}
      return pos;
    }

  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 2);
  unsigned long Mask = 0xful << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 2);
  Mask = 0xeul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Mask = 0xdul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Mask = 0xcul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 2);
  Mask = 0xbul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 2);
  Mask = 0xaul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = 0x9ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = 0x8ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = 0x7ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = 0x6ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp -  1, nbrFermionsPlus);
  Mask = 0x5ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1, nbrFermionsPlus);
  Mask = 0x4ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp, nbrFermionsPlus - 1);
  Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp, nbrFermionsPlus - 1);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp, nbrFermionsPlus);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos, nbrFermionsUp, nbrFermionsPlus);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrParticlesDownMinus = number of particles with down minus
// nbrParticlesDownPlus = number of particles with down plus
// nbrParticlesUpMinus = number of particles with up minus
// nbrParticlesUpPlus = number of particles with up plus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
								    int nbrParticlesDownMinus, int nbrParticlesDownPlus, int nbrParticlesUpMinus, int nbrParticlesUpPlus, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)
      || (nbrParticlesDownMinus < 0) || (nbrParticlesDownPlus < 0) || (nbrParticlesUpMinus < 0) || (nbrParticlesUpPlus < 0)
      || (nbrParticlesDownMinus > nbrFermions) || (nbrParticlesDownPlus > nbrFermions) || (nbrParticlesUpMinus > nbrFermions) || (nbrParticlesUpPlus > nbrFermions))
    {
      return 0l;
    }
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
      unsigned long TmpMask = ((unsigned long) (nbrParticlesDownMinus | (nbrParticlesDownPlus << 1) | (nbrParticlesUpMinus << 2) | (nbrParticlesUpPlus << 3)));
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
	      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = TmpMask << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		  (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = TmpMask << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		}
	    }
	}
      return pos;
    }

  long TmpPos;
  unsigned long Mask;
  int TmpNbrParticles;  
  for (int i = 15; i >= 0; --i)
    {
      TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3));
      TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
				    currentTotalKx + (TmpNbrParticles * currentKx),
				    currentTotalKy + (TmpNbrParticles * currentKy),
				    nbrParticlesDownMinus - (i & 1), nbrParticlesDownPlus - ((i & 2) >> 1),
				    nbrParticlesUpMinus - ((i & 4) >> 2), nbrParticlesUpPlus - ((i & 8) >> 3), pos);
      Mask = ((unsigned long) i) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
   }
  return pos;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
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
	    Count += 4l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 4l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
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
// nbrFermionsUp = current number of fermions with a spin up
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrFermionsUp)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)|| (nbrFermionsUp < 0) ||  (nbrFermionsUp > nbrFermions))
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
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrFermionsUp - 2);
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2));
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1));
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 2);
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp);
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp));
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1));
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrFermionsUp);
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
// nbrFermionsPlus = current number of fermions with a plus
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrFermionsUp, int nbrFermionsPlus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)|| (nbrFermionsUp < 0) || (nbrFermionsPlus < 0) ||  
      (nbrFermionsUp > nbrFermions) || (nbrFermionsPlus > nbrFermions))
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
	    Count += 1l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 1l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 2);
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 2);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus - 2);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus - 1);
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus -1));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp -1, nbrFermionsPlus - 2);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp, nbrFermionsPlus - 1);
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp, nbrFermionsPlus);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1, nbrFermionsPlus);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp, nbrFermionsPlus - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrFermionsUp, nbrFermionsPlus);
  return Count;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrParticlesDownMinus = number of particles with down minus
// nbrParticlesDownPlus = number of particles with down plus
// nbrParticlesUpMinus = number of particles with up minus
// nbrParticlesUpPlus = number of particles with up plus
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
										   int nbrParticlesDownMinus, int nbrParticlesDownPlus, int nbrParticlesUpMinus, int nbrParticlesUpPlus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)
      || (nbrParticlesDownMinus < 0) || (nbrParticlesDownPlus < 0) || (nbrParticlesUpMinus < 0) || (nbrParticlesUpPlus < 0)
      || (nbrParticlesDownMinus > nbrFermions) || (nbrParticlesDownPlus > nbrFermions) || (nbrParticlesUpMinus > nbrFermions) || (nbrParticlesUpPlus > nbrFermions))
    {
      return 0l;
    }
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
  
  unsigned long Tmp = 0l;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Tmp += 1l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Tmp += 1l;
	    }
	}
      return Tmp;
    }
  
  for (int i = 15; i >= 0; --i)
    {
      int TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3));
      Tmp += this->EvaluateHilbertSpaceDimension(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
						 currentTotalKx + (TmpNbrParticles * currentKx), currentTotalKy + (TmpNbrParticles * currentKy),
						 nbrParticlesDownMinus - (i & 1), nbrParticlesDownPlus - ((i & 2) >> 1),
						 nbrParticlesUpMinus - ((i & 4) >> 2), nbrParticlesUpPlus - ((i & 8) >> 3));
    }
  return Tmp;
}


// provide the default name of the file for Hilbert space storage 
//
// return value = pointer to file name (0 if no default file name exists) 

char* FermionOnSquareLatticeWithSU4SpinMomentumSpace::GetDefaultHilbertSpaceFileName()
{
  char* TmpName = new char[256];
  if ((this->SzFlag == false) && (this->PzFlag == false))
    {
      sprintf(TmpName, "fermions_klattice_su4_n_%d_nx_%d_ny_%d_kx_%d_ky_%d.hil", this->NbrFermions, this->NbrSiteX, this->NbrSiteY, this->KxMomentum, this->KyMomentum);
    }
  if ((this->SzFlag == true) && (this->PzFlag == false))
    {
      sprintf(TmpName, "fermions_klattice_su4_n_%d_nx_%d_ny_%d_kx_%d_ky_%d_sz_%d.hil", this->NbrFermions, this->NbrSiteX, this->NbrSiteY, this->KxMomentum, this->KyMomentum, this->TotalSpin);
    }
  if ((this->SzFlag == false) && (this->PzFlag == true))
    {
      sprintf(TmpName, "fermions_klattice_su4_n_%d_nx_%d_ny_%d_kx_%d_ky_%d_pz_%d.hil", this->NbrFermions, this->NbrSiteX, this->NbrSiteY, this->KxMomentum, this->KyMomentum, this->TotalIsospin);
    }
  if ((this->SzFlag == true) && (this->PzFlag == true))
    {
      sprintf(TmpName, "fermions_klattice_su4_n_%d_nx_%d_ny_%d_kx_%d_ky_%d_sz_%d_pz_%d.hil", this->NbrFermions, this->NbrSiteX, this->NbrSiteY, this->KxMomentum, this->KyMomentum, this->TotalSpin, this->TotalIsospin);
    }
  return TmpName;
}

// evaluate all the single band Hilbert spaces
//
// maxBandOccupation = maiximum occupation of a single band
// singleBandTotalKxMax = reference on the array for the maximum total Kx values 
// singleBandTotalKyMax = reference on the array for the maximum total Ky values
// singleBandHilbertDimensions = reference on the array for all the Hilbert space dimension (first index being the particle number, second index being the total Kx, third index being the total Ky)
// singleBandStates = reference on the array for all the Hilbert space basis states 

void FermionOnSquareLatticeWithSU4SpinMomentumSpace::GenerateAllSingleBandHilbertSpaces(int maxBandOccupation, int*& singleBandTotalKxMax, int*& singleBandTotalKyMax, long***& singleBandHilbertDimensions, unsigned long****& singleBandStates)
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
      TmpHoleMask |= 0x1ul << (4 * i);
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

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::EvaluateSingleBandHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy)
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

long FermionOnSquareLatticeWithSU4SpinMomentumSpace::GenerateSingleBandStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int singleBandTotalKx, int singleBandTotalKy, unsigned long* singleBandStateDescription, long pos)
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
  ULONGLONG Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) * 4);
  for (; pos < TmpPos; ++pos)
    singleBandStateDescription[pos] |= Mask;
  return this->GenerateSingleBandStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, singleBandTotalKx, singleBandTotalKy, singleBandStateDescription, pos);
};


