////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                     without Sz conservation and using Gutziller            //
//                          projection in orbital space                       //
//                                                                            //
//                        last modification : 08/12/2023                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinAllSzGutzwillerProjection.h"
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

FermionOnSphereWithSpinAllSzGutzwillerProjection::FermionOnSphereWithSpinAllSzGutzwillerProjection()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion// totalSpin = twce the total spin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSzGutzwillerProjection::FermionOnSphereWithSpinAllSzGutzwillerProjection (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
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

FermionOnSphereWithSpinAllSzGutzwillerProjection::FermionOnSphereWithSpinAllSzGutzwillerProjection(const FermionOnSphereWithSpinAllSzGutzwillerProjection& fermions)
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

FermionOnSphereWithSpinAllSzGutzwillerProjection::~FermionOnSphereWithSpinAllSzGutzwillerProjection ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinAllSzGutzwillerProjection& FermionOnSphereWithSpinAllSzGutzwillerProjection::operator = (const FermionOnSphereWithSpinAllSzGutzwillerProjection& fermions)
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

AbstractHilbertSpace* FermionOnSphereWithSpinAllSzGutzwillerProjection::Clone()
{
  return new FermionOnSphereWithSpinAllSzGutzwillerProjection(*this);
}



// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinAllSzGutzwillerProjection::GenerateStates(int nbrFermions, int posMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (posMax < (nbrFermions - 1)))
    return pos;
  
  int LzTotalMax = this->MaxTotalLz[posMax][nbrFermions];  
  
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
  TmpPos = this->GenerateStates(nbrFermions - 1, posMax - 2, totalLz - (posMax >> 1),  pos);
  Mask = 0x1ul << posMax;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, posMax - 2, totalLz - (posMax >> 1),  pos);
  Mask = 0x1ul << (posMax - 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, posMax - 2, totalLz, pos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// posMax = highest position for next particle to be placed
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereWithSpinAllSzGutzwillerProjection::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int posMax, int totalLz)
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
  long Count = 0l;
  Count += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, posMax - 2, totalLz - (posMax>>1));
  Count += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, posMax - 2, totalLz - (posMax>>1));
  Count += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, posMax - 2, totalLz);
  return Count;
}

// convert a given state from a generic basis to the gutzwiller basis
//
// state = reference on the vector to convert
// basis = pointer to the basis associated to state
// return value = converted vector

RealVector FermionOnSphereWithSpinAllSzGutzwillerProjection::GutzwillerProjection(RealVector& state, ParticleOnSphereWithSpin* basis)
{
  FermionOnSphereWithSpinAllSz* TmpSpace = (FermionOnSphereWithSpinAllSz*) basis;
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NewLzMax = 1 + (this->LzMax << 1);
      unsigned long TmpState = this->StateDescription[i];
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
      int TmpIndex = TmpSpace->FindStateIndex(TmpState, NewLzMax);
      if (TmpIndex != TmpSpace->GetHilbertSpaceDimension())
	{
	  TmpVector[i] = state[TmpIndex];
	}
    }
  return TmpVector;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix FermionOnSphereWithSpinAllSzGutzwillerProjection::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
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
 	  FermionOnSphereWithSpinAllSzGutzwillerProjection TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize);
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
      FermionOnSphereWithSpinAllSzGutzwillerProjection TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1, true);
      int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	    TmpEntanglementMatrix.AddToMatrixElement(i, 0, groundState[MinIndex + i]);
	}
      return TmpEntanglementMatrix;
    }

  FermionOnSphereWithSpinAllSzGutzwillerProjection TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
 
  FermionOnSphereWithSpinAllSzGutzwillerProjection TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize);
 
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

RealMatrix FermionOnSphereWithSpinAllSzGutzwillerProjection::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, 
														 bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  int nbrOrbitalA = this->LzMax + 1;
  int nbrOrbitalB = this->LzMax + 1;  

  if (nbrParticleSector == 0)
    {
      if (lzSector == 0)
        {
          FermionOnSphereWithSpinAllSzGutzwillerProjection TmpHilbertSpace(this->NbrFermions, this->TotalLz - lzSector, this->LzMax);
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
          FermionOnSphereWithSpinAllSzGutzwillerProjection TmpDestinationHilbertSpace(nbrParticleSector, lzSector, this->LzMax);
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

  FermionOnSphereWithSpinAllSzGutzwillerProjection SubsytemSpace(nbrParticleSector, lzSector, this->LzMax);
  FermionOnSphereWithSpinAllSzGutzwillerProjection ComplementarySubsytemSpace(ComplementaryNbrParticles, this->TotalLz - lzSector, this->LzMax);

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
