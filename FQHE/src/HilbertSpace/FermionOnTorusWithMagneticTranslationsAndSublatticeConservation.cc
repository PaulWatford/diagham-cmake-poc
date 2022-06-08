////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//           and conservation of the number of particles per sublattice       //
//                                                                            //
//                        last modification : 05/06/2022                      //
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
#include "HilbertSpace/FermionOnTorusWithMagneticTranslationsAndSublatticeConservation.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include <math.h>
#include <limits.h>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// default constructor
// 

FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::FermionOnTorusWithMagneticTranslationsAndSublatticeConservation ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
// evenMomentumNbrParticles = number of particles occupying orbitals with an even momentum

FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::FermionOnTorusWithMagneticTranslationsAndSublatticeConservation (int nbrFermions, int maxMomentum, 
																  int xMomentum, int yMomentum, int evenMomentumNbrParticles)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;  
  this->EvenMomentumNbrParticles = evenMomentumNbrParticles;
  
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = 0x1ul;
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= 0x1ul;
    }
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->MaxMomentum + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->MaxMomentum + 1) >> 1);
#endif
  if ((this->MaxMomentum & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  long TmpDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum);
  cout << "Max dimension: " << TmpDimension << endl;
  if (TmpDimension>INT_MAX)
    cout << "Max-Dimension surpasses integer representation..."<<endl;
  this->StateDescription = new unsigned long [TmpDimension];
  this->StateMaxMomentum = new int [TmpDimension];
  this->LargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->MaxMomentum - 1, 0l, 0, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Actual dimension: " << this->LargeHilbertSpaceDimension << endl;

  this->Flag.Initialize();
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
      int UsedMemory = 0;
      UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
      UsedMemory += this->NbrMomentum * sizeof(int);
      //  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
      UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
  else
    {
      this->LookUpTableShift = 0;
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::FermionOnTorusWithMagneticTranslationsAndSublatticeConservation(const FermionOnTorusWithMagneticTranslationsAndSublatticeConservation& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->EvenMomentumNbrParticles = fermions.EvenMomentumNbrParticles;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;
}

// destructor
//

FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::~FermionOnTorusWithMagneticTranslationsAndSublatticeConservation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithMagneticTranslationsAndSublatticeConservation& FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::operator = (const FermionOnTorusWithMagneticTranslationsAndSublatticeConservation& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->EvenMomentumNbrParticles = fermions.EvenMomentumNbrParticles;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::Clone()
{
  return new FermionOnTorusWithMagneticTranslationsAndSublatticeConservation(*this);
}

// generate all states corresponding to the constraints (without taking into the canonical form) 
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentYMomentum = current value of the momentum in the y direction
// currentEvenMomentumNbrParticles = current number of particles occupying orbitals with an even momentum
// return value = position from which new states have to be stored

long FermionOnTorusWithMagneticTranslationsAndSublatticeConservation::RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, long pos, int currentYMomentum, int currentEvenMomentumNbrParticles)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)))
    return pos;
  if (nbrFermions == 1)
    {
      int i = this->YMomentum - (currentYMomentum % this->MaxMomentum);
      if (i < 0)
 	i += this->MaxMomentum;
      for (; i <= currentMaxMomentum; i += this->MaxMomentum)
	{
	  if ((((i + 1) & 1) + currentEvenMomentumNbrParticles) == this->EvenMomentumNbrParticles)
	    {
	      this->StateDescription[pos] = 0x1ul << i;
	      this->StateMaxMomentum[pos] = maxMomentum;
	      ++pos;
	    }
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  long TmpPos = this->RawGenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, currentYMomentum + currentMaxMomentum, currentEvenMomentumNbrParticles + ((currentMaxMomentum + 1) & 1));
  unsigned long Mask = 0x1ul << currentMaxMomentum;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentYMomentum, currentEvenMomentumNbrParticles);
  else
    return this->RawGenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentYMomentum, currentEvenMomentumNbrParticles);
}

