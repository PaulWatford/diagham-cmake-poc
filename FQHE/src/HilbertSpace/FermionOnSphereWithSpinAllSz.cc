////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                            sign precalculation table                       //
//                                                                            //
//                        last modification : 12/12/2005                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/FactorialCoefficient.h"

#include <cmath>
#include <bitset>
#include <cstdlib>

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

FermionOnSphereWithSpinAllSz::FermionOnSphereWithSpinAllSz()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion// totalSpin = twce the total spin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSz::FermionOnSphereWithSpinAllSz (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TargetSpace = this;

#ifdef  __64_BITS__
  if (this->NbrLzValue>32)
    {
      cout<<"Cannot represent the system size requested in a single word: only LzMax<=31 supported"<<endl;
      exit(1);
    }
#else
  if (this->NbrLzValue>16)
    {
      cout<<"Cannot represent the system size requested in a single word: only LzMax<=15 supported"<<endl;
      exit(1);
    }
#endif

  // temporary space to accelerate state generation
  this->MaxTotalLz = new int*[2*NbrLzValue];
  for (int i=0; i<2*NbrLzValue; ++i)
    {
      MaxTotalLz[i] = new int[NbrFermions+1];
      for (int f=0; f<=NbrFermions; ++f)
	{
	  MaxTotalLz[i][f] = 0;
	  for (int n=0; (n<f); ++n)	  
	    MaxTotalLz[i][f] += (i-n)>>1;
 	}
    }
    
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, (this->LzMax<<1)+1, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
  this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, (this->LzMax<<1)+1, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0x0l);
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  // clean up temporary space
  for (int i=0; i<2*NbrLzValue; ++i)
    delete [] MaxTotalLz[i];
  delete [] MaxTotalLz;

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      PrintMemorySize(cout, UsedMemory)<<endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      PrintMemorySize(cout,UsedMemory)<<endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSpinAllSz::FermionOnSphereWithSpinAllSz(const FermionOnSphereWithSpinAllSz& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->TargetSpace = this;
}

// destructor
//

FermionOnSphereWithSpinAllSz::~FermionOnSphereWithSpinAllSz ()
{
/*
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < (2 * this->NbrLzValue); ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
*/
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinAllSz& FermionOnSphereWithSpinAllSz::operator = (const FermionOnSphereWithSpinAllSz& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->TargetSpace = this;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinAllSz::Clone()
{
  return new FermionOnSphereWithSpinAllSz(*this);
}


// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinAllSz::AduAd (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];  

  m <<= 1;
  ++m;
  n <<= 1;

  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  // perform annihilation operators
  TmpState &= ~(((unsigned long) (0x1)) << n);
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  unsigned long TmpState2=TmpState;
  if (NewHighestBit == n)
    while ((NewHighestBit > 0)&&((TmpState >> NewHighestBit) == 0))
      --NewHighestBit;
  // create particle at m
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewHighestBit)
    {
      NewHighestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m);
//  return this->FindStateIndex(TmpState, NewHighestBit);
  return this->FindStateIndex(TmpState,NewHighestBit);
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinAllSz::AddAu (int index, int m, int n, double& coefficient)
{
    int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];  

  m <<= 1;
  n <<= 1;
  ++n;

  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  // perform annihilation operators
  TmpState &= ~(((unsigned long) (0x1)) << n);
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  if (NewHighestBit == n)
    while ((NewHighestBit > 0)&&((TmpState >> NewHighestBit) == 0))
      --NewHighestBit;
  // create particle at m
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewHighestBit)
    {
      NewHighestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m);
//  return this->FindStateIndex(TmpState, NewHighestBit);
  return this->FindStateIndex(TmpState,NewHighestBit);
}
  

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinAllSz::GenerateStates(int nbrFermions, int posMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (posMax < (nbrFermions - 1)))
    return pos;
  
//   if ((nbrFermions == 0) && (totalLz == 0)) // this could occur only if we assigned two fermions at once
//     {
//       cout << "0 fermions"<<endl;
//       this->StateDescription[pos] = 0x0ul;
//       this->StateHighestBit[pos++] = 0;
//       return pos;
//     }

  int LzTotalMax = this->MaxTotalLz[posMax][nbrFermions];  
//   int LzTotalMax = 0;
//   for (int n=0; n<nbrFermions; ++n) LzTotalMax += (currentPosMax-n)/2;
  
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
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// posMax = highest position for next particle to be placed
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereWithSpinAllSz::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int posMax, int totalLz)
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


// Project the state from the tunneling space (all Sz's)
// to the space with the fixed projection of Sz (given by SzValue)
//
// state = state that needs to be projected
// su2Space = the subspace onto which the projection is carried out
// SzValue = the desired value of Sz

RealVector FermionOnSphereWithSpinAllSz::ForgeSU2FromTunneling(RealVector& state, FermionOnSphereWithSpin& su2Space, int SzValue)
{
  RealVector FinalState(su2Space.GetHilbertSpaceDimension(), true);
  int counter=0;
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
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
       FinalState[su2Space.FindStateIndex(this->StateDescription[j], this->StateHighestBit[j])] += state[j];
	}
    }
  cout<<"Nbr of stored components = "<<counter<<endl;
  //FinalState /= FinalState.Norm();
  return FinalState;  
}


// Project the state from the tunneling space (all Sz's)
// to the U(1) space (u1Space)
//
// state = state that needs to be projected
// u1Space = the subspace onto which the projection is carried out

RealVector  FermionOnSphereWithSpinAllSz::ForgeU1FromTunneling(RealVector& state, FermionOnSphere& u1Space)
{
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  int counter=0;
  int rejected=0;

  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      int TmpPos = this->LzMax << 1;
      int TmpUpValue, TmpDownValue, TmpUpDownValue;
      int nbr_X=0;
      while (TmpPos >= 0)
	{
 	  TmpUpValue = (TmpState>>(TmpPos+1)) & 0x1ul ;
 	  TmpDownValue = (TmpState>>TmpPos) & 0x1ul ;
	  if ((TmpUpValue==1)&&(TmpDownValue==1) )
		++nbr_X;
	  TmpPos -= 2;
	}
      if  (nbr_X==0) //If there is no orbital filled with both Up and Down, proceed...
	{ 
 	  ++counter;
	  TmpPos = this->LzMax << 1;
	  unsigned long TmpState2 = 0x0ul;
	  int tmp_counter=this->LzMax;
	  while (TmpPos >= 0)
	  {
	      TmpUpDownValue = ((TmpState >> TmpPos) & ((unsigned long) 0x3));
	      if (TmpUpDownValue == 0x1l)
 		TmpState2 |=  ((unsigned long)0x1 << tmp_counter); //"d ";
	      else if (TmpUpDownValue == 0x2l)
 		TmpState2 |=  ((unsigned long)0x1 << tmp_counter); //"u ";
	      else if (TmpUpDownValue == 0x3l)
 		TmpState2 |=  ((unsigned long)0x1 << tmp_counter); //"X ";
	      else  TmpState2 |=  ((unsigned long)0x0 << tmp_counter); //"0 ";
		TmpPos -= 2;
		--tmp_counter;
	  }
          TmpPos=this->LzMax;
	  while ((TmpState2 >> TmpPos) == 0x0ul)
		--TmpPos;

	  if (u1Space.FindStateIndex( TmpState2, TmpPos)==0)
	   {	
	     if (TmpState2==u1Space.StateDescription[0])
		FinalState[0] += state[j]*pow(1.0/sqrt(2.0), this->NbrFermions);
	     else 
		rejected++;

	   }
	  else
	   FinalState[u1Space.FindStateIndex( TmpState2, TmpPos)]+=state[j]*pow(1.0/sqrt(2.0), this->NbrFermions);
	  
	} //End if (nbr_X==0)

    } //End loop over HilbertSpace

  cout<<"Nbr of kept components = "<<counter<<endl;
  cout<<"Rejected "<<rejected<<endl; 
  //FinalState /= FinalState.Norm();
  return FinalState;  
}

// Calculate mean value <Sx> in a given state
//
// state = given state

double FermionOnSphereWithSpinAllSz::MeanSxValue(RealVector& state)
{
  double Coefficient;
  int Dim = this->HilbertSpaceDimension;

  RealVector FinalState(Dim, true);

  for (int i = 0; i < Dim; ++i)    
   for (int j=0; j<=this->LzMax; ++j)
    {
      int Index=this->AduAd(i,j,j,Coefficient);
      if ( (Index<Dim) && (Coefficient != 0.0))
	FinalState[Index]+=state[i]*Coefficient*0.5;
      Index=this->AddAu(i,j,j,Coefficient);
      if ( (Index<Dim) && (Coefficient != 0.0))
	FinalState[Index]+=state[i]*Coefficient*0.5;
     }

  return (FinalState*state);
}

// Calculate mean value <Sz> in a given state
//
// state = given state

double FermionOnSphereWithSpinAllSz::MeanSzValue(RealVector& state)
{
  double Coefficient;
  int Dim = this->HilbertSpaceDimension;

  RealVector FinalState(Dim, true);

  for (int i = 0; i < Dim; ++i)    
      for (int j=0; j<=this->LzMax; ++j)
	{
	 FinalState[i]+=state[i]*0.5*this->AduAu(i,j);
	 FinalState[i]-=state[i]*0.5*this->AddAd(i,j);
	}

  return (FinalState*state);
}

// Artificially extend a state of a U(1) Hilbert space to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// u1space = U(1) space of the input state
// return value = input state expression in the SU(2) basis

ComplexVector FermionOnSphereWithSpinAllSz::U1ToSU2AllSz(ComplexVector& state, FermionOnSphere& u1space)
{
  ComplexVector FinalState(this->GetHilbertSpaceDimension(), true);
  FinalState.ClearVector();
  int U1Dimension=state.GetVectorDimension();
  unsigned long U1BasisState;
  unsigned long SU2BasisState;
  unsigned long TmpState;

  for (int j = 0; j < U1Dimension; ++j)    
    {
      U1BasisState = u1space.StateDescription[j];
      SU2BasisState = 0x0ul;

      // Print U1BasisState
//      cout << "U1BasisState number " << j << ": " << bitset<sizeof(unsigned long)*8>(U1BasisState) << " -> SU2BasisState: ";
      for(int k = 0; k < this->NbrLzValue; ++k)
      {
	TmpState = U1BasisState % 2; 
	SU2BasisState += TmpState << (2*k+1);
	U1BasisState >>= 1;
      }
      // Print SU2BasisState
//      cout << bitset<sizeof(unsigned long)*8>(SU2BasisState) << endl;
      
      int TmpLzMax = this->NbrLzValue << 1;
      while ((SU2BasisState >> TmpLzMax) == 0x0ul)
	--TmpLzMax; //TmpLzMax = maximum Lz with non zero occupation in TmpState

//      int index = this->CarefulFindStateIndex(SU2BasisState,-1);
      int index = this->FindStateIndex(SU2BasisState,TmpLzMax);
      FinalState[index]=state[j];
    }
  
//  FinalState /= FinalState.Norm();
  return FinalState;  
}

// Artificially extend a state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// su2space = SU(2) space with fixed sz of the input state
// return value = input state expression in the SU(2) basis

ComplexVector FermionOnSphereWithSpinAllSz::SU2ToSU2AllSz(ComplexVector& state, FermionOnSphereWithSpin& su2space)
{
  ComplexVector FinalState(this->GetHilbertSpaceDimension(), true);
  FinalState.ClearVector();
  int SU2Dimension=state.GetVectorDimension();
  unsigned long SU2BasisState;
  unsigned long TmpState;

  for (int j = 0; j < SU2Dimension; ++j)    
    {
      SU2BasisState = su2space.StateDescription[j];
     
      int TmpLzMax = this->NbrLzValue << 1;
      while ((SU2BasisState >> TmpLzMax) == 0x0ul)
	--TmpLzMax; //TmpLzMax = maximum Lz with non zero occupation in TmpState

//      int index = this->CarefulFindStateIndex(SU2BasisState,-1);
      int index = this->FindStateIndex(SU2BasisState,TmpLzMax);
      FinalState[index]=state[j];
    }
  
//  FinalState /= FinalState.Norm();
  return FinalState;  
}

// convert a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void FermionOnSphereWithSpinAllSz::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, 
							 long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrFermions];
  int* TmpSpinIndices = new int [this->NbrFermions];
  int* TmpSpinIndices2 = new int [this->NbrFermions];
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      unsigned long Tmp;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  Tmp = (TmpState >> (j << 1)) & 0x3ul;;
	  if ((Tmp & 0x2ul) != 0x0ul) // If there is a particle with momentum j and spin down in the ith state of the basis 
	    {
	      TmpMomentumIndices[TmpIndex] = j; // An array which gathers all momenta
	      TmpSpinIndices[TmpIndex] = 1; // An array which gathers all spins: 1 = down
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x1ul) != 0x0ul) // If there is a particle with momentum j and spin up in the ith state of the basis  
	    {
	      TmpMomentumIndices[TmpIndex] = j;// An array which gathers all momenta
	      TmpSpinIndices[TmpIndex] = 0;// An array which gathers all spins: 0 = up
	      ++TmpIndex;
	    }	  
	}	
	// Print state an arrays
//	cout << "Basis state number = " <<  i << ": " << bitset<sizeof(unsigned long)*8>(this->StateDescription[i]) << endl;
//	cout << "Momentum indices: " << endl;
//	for( int k = 0; k < this->NbrFermions ; ++k) 
//	  cout << TmpMomentumIndices[k] << " ";
//	cout <<endl;
//	cout << "SpinIndices: " << endl;
//	for( int k = 0; k < this->NbrFermions ; ++k) 
//	  cout << TmpSpinIndices[k] << " ";
//	cout <<endl;

	// initialState[i]: coefficient of the state to be transformed on the ith vector of the basis, has to multiply the final state
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSpinIndices, TmpSpinIndices2, oneBodyBasis);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSpinIndices;
  delete[] TmpSpinIndices2;
}

// recursive part of the convertion from a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSpinIndices = array that gives the spin dressing the initial n-body state
// currentSpinIndices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector

void FermionOnSphereWithSpinAllSz::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
								int position, int* momentumIndices, int* initialSpinIndices, int* currentSpinIndices, ComplexMatrix* oneBodyBasis) 
{
//   cout << position << " : " << endl;
//   for (int i = 0; i < position; ++i)
//     cout << currentSpinIndices[i] << " ";
//   cout << endl;
  if (position == this->NbrFermions)
    {
      unsigned long TmpState = 0x0ul;
      unsigned long TmpState2;
      unsigned long Mask = 0x0ul;
      unsigned long MaskSign = 0x0ul;
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  Mask = 0x1ul << ((momentumIndices[i] << 1) + currentSpinIndices[i]); // Mask = 00...0100...0 : one fermion state in the second quantized basis
	  if ((TmpState & Mask) != 0x0ul)
	    return;
	  // SignMask computation -----------------------------------
	  TmpState2 = TmpState & (Mask - 0x1ul); 
#ifdef __64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= (TmpState2 >> 16);
	  TmpState2 ^= (TmpState2 >> 8);
	  TmpState2 ^= (TmpState2 >> 4);
	  TmpState2 ^= (TmpState2 >> 2);
	  MaskSign ^= (TmpState2 ^ (TmpState2 >> 1)) & 0x1ul;
	  // End of SignMask computation -----------------------------------

	  TmpState |= Mask; //set bit corresponding to the current fermion state to 1 in TmpState
	}
      int TmpLzMax = this->NbrLzValue << 1;
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax; //TmpLzMax = maximum Lz with non zero occupation in TmpState
      int Index = this->FindStateIndex(TmpState, TmpLzMax);
      if (Index < this->HilbertSpaceDimension)
	{
	  if (MaskSign == 0ul)
	    {
	      targetState[Index] += coefficient;
	    }
	  else
	    {
	      targetState[Index] -= coefficient;
	    }
	}
      return;      
    }
  else
    {
      currentSpinIndices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSpinIndices[position]][1]), position + 1, momentumIndices, initialSpinIndices, currentSpinIndices, oneBodyBasis);
      currentSpinIndices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSpinIndices[position]][0]), position + 1, momentumIndices, initialSpinIndices, currentSpinIndices, oneBodyBasis);
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix FermionOnSphereWithSpinAllSz::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
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
	  double TmpValue = 0;
 	  FermionOnSphereWithSpinAllSz TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize);
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
	      int TmpLzMax = 2 * this->LzMax + 1;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
                {
                   TmpNbrNonZeroElements++;
		   TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex, groundState[TmpPos]);	
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

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = 2 * this->LzMax + 1;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
              TmpNbrNonZeroElements++;
              TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, groundState[TmpPos]);
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

RealMatrix FermionOnSphereWithSpinAllSz::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, 
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
      if (lzSector == this->TotalLz)
        {
          FermionOnSphereWithSpinAllSz TmpDestinationHilbertSpace(nbrParticleSector, lzSector, this->LzMax);
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

  FermionOnSphereWithSpinAllSz SubsytemSpace(nbrParticleSector, lzSector, this->LzMax);
  FermionOnSphereWithSpinAllSz ComplementarySubsytemSpace(ComplementaryNbrParticles, this->TotalLz - lzSector, this->LzMax);

  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);

  long TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore(0, ComplementarySubsytemSpace.GetHilbertSpaceDimension(),
											    &ComplementarySubsytemSpace, &SubsytemSpace, 
											    groundState, &TmpEntanglementMatrix, removeBinomialCoefficient);
  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
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

long FermionOnSphereWithSpinAllSz::EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, 
											   ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealMatrix* entanglementMatrix, 
											   bool removeBinomialCoefficient)
{
  FermionOnSphereWithSpinAllSz* SubsytemSpace = (FermionOnSphereWithSpinAllSz*) ((FermionOnSphereWithSpinAllSz*) destinationHilbertSpace)->Clone();
  FermionOnSphereWithSpinAllSz* ComplementarySubsytemSpace = (FermionOnSphereWithSpinAllSz*) ((FermionOnSphereWithSpinAllSz*) complementaryHilbertSpace)->Clone();

  long TmpNbrNonZeroElements = 0;
  double TmpInvBinomial = 1.0;
  if (removeBinomialCoefficient == false)
    {
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      TmpInvBinomial = sqrt(1.0 / (TmpBinomial( this->NbrFermions , SubsytemSpace->NbrFermions))); 
    }

  int MaxIndex = minIndex + nbrIndex;
  for (int MinIndex = minIndex; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = ComplementarySubsytemSpace->StateDescription[MinIndex];
      for (int j = 0; j < SubsytemSpace->HilbertSpaceDimension; ++j)
        {
          unsigned long TmpState2 = SubsytemSpace->StateDescription[j];
          if ((TmpState & TmpState2) == 0x0ul)
            {
              int TmpLzMax = (this->LzMax << 1) + 1;
              unsigned long TmpState3 = TmpState | TmpState2;
              while ((TmpState3 >> TmpLzMax) == 0x0ul)
                --TmpLzMax;
              int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
              if (TmpPos != this->HilbertSpaceDimension)
                {
                  double Coefficient = TmpInvBinomial;
                  unsigned long Sign = 0x0ul;
                  int Pos2 = (SubsytemSpace->LzMax << 1) + 1;
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
                  entanglementMatrix->SetMatrixElement(j, MinIndex, Coefficient*groundState[TmpPos]);
                }
            }
        }
    }
  
  delete SubsytemSpace;
  delete ComplementarySubsytemSpace;
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

RealMatrix& FermionOnSphereWithSpinAllSz::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
															     int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
															     int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealMatrix& entanglementMatrix)
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

