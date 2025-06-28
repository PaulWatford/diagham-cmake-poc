////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on square lattice with spin            //
//         in momentum space using a simplified version more in line with     //
//                         SU(3), SU(6) or SU(12) codes                       //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpace.h"
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

#include <math.h>
#include <cstdlib>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
//

FermionOnSquareLatticeWithSU2SpinMomentumSpace::FermionOnSquareLatticeWithSU2SpinMomentumSpace()
{
}
  
// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU2SpinMomentumSpace::FermionOnSquareLatticeWithSU2SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
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
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  timeval TotalEndingTime;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	       ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));                   
  cout << "Hilbert space dimension computed in " << Dt << "s" << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      gettimeofday (&(TotalStartingTime), 0);
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      gettimeofday (&(TotalEndingTime), 0);
      Dt = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));                   
      cout << "Hilbert space generated in " << Dt << "s" << endl;
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
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU2SpinMomentumSpace::FermionOnSquareLatticeWithSU2SpinMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->NbrFermionsUp = nbrSpinUp;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->TotalSpin = (this->NbrFermionsUp - this->NbrFermionsDown);
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);
  if (this->NbrFermions <= (this->NbrSiteX * this->NbrSiteY))
    {
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp);
    }
  else
    {
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimensionHoles((2 * this->NbrSiteX * this->NbrSiteY) - this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, (this->NbrSiteX * (this->NbrSiteX - 1)) * this->NbrSiteY, (this->NbrSiteY * (this->NbrSiteY - 1)) * this->NbrSiteX, (this->NbrSiteX * this->NbrSiteY) - this->NbrFermionsUp);      
    }
  timeval TotalEndingTime;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	       ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));                   
  cout << "Hilbert space dimension computed in " << Dt << "s" << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = 0l;
      gettimeofday (&(TotalStartingTime), 0);
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
      gettimeofday (&(TotalEndingTime), 0);
      Dt = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));                   
      cout << "Hilbert space generated in " << Dt << "s" << endl;
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

FermionOnSquareLatticeWithSU2SpinMomentumSpace::FermionOnSquareLatticeWithSU2SpinMomentumSpace(const FermionOnSquareLatticeWithSU2SpinMomentumSpace& fermions)
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

FermionOnSquareLatticeWithSU2SpinMomentumSpace::~FermionOnSquareLatticeWithSU2SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU2SpinMomentumSpace& FermionOnSquareLatticeWithSU2SpinMomentumSpace::operator = (const FermionOnSquareLatticeWithSU2SpinMomentumSpace& fermions)
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU2SpinMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU2SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeWithSU2SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 1));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & 0x2l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",+)";
      if ((Tmp & 0x1l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",-)";
    }
  Str << "]";
  return Str;
}

// evaluate Hilbert space dimension with a fixed number of bosons with spin up
//
// initialState = state to be projected
// finalSpace = space in which the projected state is
// finalState = state after projection
// bandIndex = index of the project to project onto

void FermionOnSquareLatticeWithSU2SpinMomentumSpace::ProjectOntoSingleBand(ComplexVector* initialState, ParticleOnSphere* finalSpace, ComplexVector* finalState, bool bandIndex)
{
  FermionOnSquareLatticeMomentumSpace* TmpSpace = (FermionOnSquareLatticeMomentumSpace*) finalSpace;
  for (long i = 0l; i < TmpSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpInputState = TmpSpace->StateDescription[i];
      unsigned long TmpOutputState = 0x0ul;
      for (int j = 0; j < 64; ++j)
	{
	  TmpOutputState |= ((TmpInputState >> j) & 0x1ul) << (j << 1);
	}
      TmpOutputState <<= bandIndex;
#ifdef __64_BITS__      
      int TmpLzMax = 63;
#else
      int TmpLzMax = 31;
#endif
      while ((TmpLzMax > 0) && (((TmpOutputState >> TmpLzMax) & 0x1ul) == 0x0ul))
	{
	  TmpLzMax--;	  
	}
      int TmpIndex = this->FindStateIndex(TmpOutputState, TmpLzMax);
      if (TmpIndex < this->HilbertSpaceDimension)
	{
	  (*finalState)[i] = (*initialState)[TmpIndex];
	}
      else
	{
	  (*finalState)[i] = 0.0;
	}
    }
}

