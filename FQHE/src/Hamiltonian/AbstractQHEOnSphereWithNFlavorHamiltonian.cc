////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of generic N-flavor Hamiltonian for particles on a           //
//                                    sphere                                  //
//                                                                            //
//                           class author: Sahana Das                         //
//                                                                            //
//                        last modification : 27/05/2026                      //
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
#include "Hamiltonian/AbstractQHEOnSphereWithNFlavorHamiltonian.h"
#include "HilbertSpace/FermionOnSphereWithNFlavor.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"


#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;



AbstractQHEOnSphereWithNFlavorHamiltonian::AbstractQHEOnSphereWithNFlavorHamiltonian()
{
  this->NbrFlavors = 0;
  this->NbrChannels = 0;

  this->PseudoPotentials = 0;
  this->OneBodyPotentials = 0;

  this->InteractionFactorsIntra = 0;
  this->InteractionFactorsInter = 0;

  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;

  this->NbrIntraSectorIndicesPerSum = 0;
  this->IntraSectorIndicesPerSum = 0;

  this->NbrInterSectorIndicesPerSum = 0;
  this->InterSectorIndicesPerSum = 0;

  this->OneBodyInteractionFactors = 0;
  this->OneBodyInteractionFactorsIntra = 0;
  this->OneBodyInteractionFactorsInter = 0;

  this->M1IntraValue = 0;
  this->M2IntraValue = 0;
  this->NbrM12IntraIndices = 0;
  this->M3IntraValues = 0;
  this->NbrM3IntraValues = 0;

  this->M1InterValue = 0;
  this->M2InterValue = 0;
  this->NbrM12InterIndices = 0;
  this->M3InterValues = 0;
  this->NbrM3InterValues = 0;

  this->FastMultiplicationStep = 0;
  this->FastMultiplicationFlag = false;
}


// virtual destructor
//

AbstractQHEOnSphereWithNFlavorHamiltonian::~AbstractQHEOnSphereWithNFlavorHamiltonian()
{
  // ---- One-body ----
  if (this->OneBodyInteractionFactors != 0)
  {
    for (int a = 0; a < this->NbrFlavors; ++a)
    {
      if (this->OneBodyInteractionFactors[a] != 0)
      {
        for (int b = 0; b < this->NbrFlavors; ++b)
        {
          // deleting null is safe, but guard the row pointer
          if (this->OneBodyInteractionFactors[a][b] != 0)
            delete[] this->OneBodyInteractionFactors[a][b];
        }
        delete[] this->OneBodyInteractionFactors[a];
      }
    }
    delete[] this->OneBodyInteractionFactors;
  }

  // ---- Intra-channel two-body ----
  if (this->InteractionFactorsIntra != 0)
  {
    int NChannels = this->NbrFlavors * (this->NbrFlavors + 1) / 2;
    for (int c = 0; c < NChannels; ++c)
    {
      if (this->InteractionFactorsIntra[c] != 0)
      {
        for (int j = 0; j < this->NbrIntraSectorSums; ++j)
        {
          if (this->InteractionFactorsIntra[c][j] != 0)
            delete[] this->InteractionFactorsIntra[c][j];
        }
        delete[] this->InteractionFactorsIntra[c];
      }
    }
    delete[] this->InteractionFactorsIntra;
  }

  // ---- Inter-channel two-body ----
  if (this->InteractionFactorsInter != 0)
  {
    int NChannels = this->NbrFlavors * (this->NbrFlavors + 1) / 2;
    for (int c = 0; c < NChannels; ++c)
    {
      if (this->InteractionFactorsInter[c] != 0)
      {
        for (int j = 0; j < this->NbrInterSectorSums; ++j)
        {
          if (this->InteractionFactorsInter[c][j] != 0)
            delete[] this->InteractionFactorsInter[c][j];
        }
        delete[] this->InteractionFactorsInter[c];
      }
    }
    delete[] this->InteractionFactorsInter;
  }
}



bool AbstractQHEOnSphereWithNFlavorHamiltonian::IsHermitian()
{
  return false;
}

bool AbstractQHEOnSphereWithNFlavorHamiltonian::IsConjugate()
{
  return false;
}


RealVector& AbstractQHEOnSphereWithNFlavorHamiltonian::
LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination,
                    int firstComponent, int nbrComponent)
{
  //std::cout << "NbrFlavors = " << this->NbrFlavors << std::endl;
  //cout << "Entering AbstractQHEOnSphereWithNFlavorHamiltonian::LowLevelAddMultiply()  "  << endl;
  //cout << "source = " << vSource << " Dest: " << vDestination  << endl;
  //cout << "firstComponent = " << firstComponent << " :: nbrComponent = " << nbrComponent << endl;

  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;

  if (this->FastMultiplicationFlag == false)
  {
    int Index;
    int* TmpIndices;
    double* TmpInteractionFactor;
    double Coefficient3;

    ParticleOnSphereWithNFlavor* TmpParticles = (ParticleOnSphereWithNFlavor*) this->Particles->Clone();

    for (int i = firstComponent; i < LastComponent; ++i)
    {
      /* ================================
         INTRA SECTOR
         ================================ */
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
      {
      
        int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
        TmpIndices = this->IntraSectorIndicesPerSum[j];
        //cout << " Intra sector sum " << j << " Lim = " << Lim << endl;

        for (int i1 = 0; i1 < Lim; i1 += 2)
        {
          for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
          {
            int channel = Channel(sigma, sigma);
            Coefficient3 = TmpParticles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], sigma, sigma);

            if (Coefficient3 != 0.0)
            {
              TmpInteractionFactor = &(this->InteractionFactorsIntra[channel][j][(i1 * Lim) >> 2]);
              //cout << " Intra sector sum " << j << ", Lim = " << Lim << ", i1 " << i1 << ", Coefficient " << sigma << " = " << Coefficient3 << " Interaction factor = " << *TmpInteractionFactor << endl;

              Coefficient3 *= vSource[i];

              for (int i2 = 0; i2 < Lim; i2 += 2)
              {
                Index = TmpParticles->AdsigmaAdsigma( TmpIndices[i2], TmpIndices[i2 + 1], sigma, sigma, Coefficient);

                if (Index < Dim)
                  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;

                ++TmpInteractionFactor;
              }
            }
          }
        }
      }
      //cout << " Intra sector vDestination = " << vDestination << endl;

      /* ================================
         INTER SECTOR
         ================================ */

      for (int j = 0; j < this->NbrInterSectorSums; ++j)
      {
        int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
        TmpIndices = this->InterSectorIndicesPerSum[j];

        for (int i1 = 0; i1 < Lim; i1 += 2)
        {
          for (int s1 = 0; s1 < NbrFlavors; ++s1)
          for (int s2 = s1 + 1; s2 < NbrFlavors; ++s2)
          {
            int channel = Channel(s1, s2);

            Coefficient3 =
              TmpParticles->AsigmaAsigma(i,
                                         TmpIndices[i1],
                                         TmpIndices[i1 + 1],
                                         s1, s2);

            if (Coefficient3 != 0.0)
            {
              TmpInteractionFactor =
                &(this->InteractionFactorsInter[channel][j]
                  [(i1 * Lim) >> 2]);

              Coefficient3 *= vSource[i];

              for (int i2 = 0; i2 < Lim; i2 += 2)
              {
                Index =
                  TmpParticles->AdsigmaAdsigma(
                    TmpIndices[i2],
                    TmpIndices[i2 + 1],
                    s1, s2,
                    Coefficient);

                if (Index < Dim)
                  vDestination[Index] +=
                    Coefficient *
                    (*TmpInteractionFactor) *
                    Coefficient3;

                ++TmpInteractionFactor;
              }
            }
          }
        }
      }   // inter sector sums
      //cout << " Inter sector vDestination = " << vDestination << endl;

    }

    /* ================================
       ONE BODY
       ================================ */

    if (this->OneBodyInteractionFactors != 0)
    {
      double TmpDiagonal = 0.0;
      for (int i = firstComponent; i < LastComponent; ++i)
      {
        TmpDiagonal = 0.0;

      /* -------- diagonal ---------- */

        for (int sigma = 0; sigma < NbrFlavors; ++sigma)
          for (int j = 0; j <= this->LzMax; ++j)
            TmpDiagonal += this->OneBodyInteractionFactors[sigma][sigma][j] * TmpParticles->AdsigmaAsigma(i, j, sigma);

        vDestination[i] += (this->HamiltonianShift + TmpDiagonal) * vSource[i];
        //cout << "One body " << i << " = " << TmpDiagonal << endl;

      /* -------- off-diagonal tunneling ---------- */

        int Index;
        double Coefficient;
        int Dim = this->Particles->GetHilbertSpaceDimension();
    
        for (int m = 0; m <= this->LzMax; ++m)
          for (int a = 0; a < NbrFlavors; ++a)
            for (int b = 0; b < NbrFlavors; ++b)
              if ((a != b) && (this->OneBodyInteractionFactors[a][b] != 0))
		{
		  Index = TmpParticles->AdsigmaAsigma(i, m, m, a, b, Coefficient);
		  
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->OneBodyInteractionFactors[a][b][m] * vSource[i];
		}
      }
      }
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}

	  //cout << " Final vDestination = " << vDestination << " Shift = " << this->HamiltonianShift << endl;

    delete TmpParticles;
  }

  else
  {
    if (this->FastMultiplicationStep == 1)
    {
      int* TmpIndexArray;
      double* TmpCoefficientArray; 
      int j;
      int TmpNbrInteraction;
      int k = firstComponent;
  
      firstComponent -= this->PrecalculationShift;
      LastComponent  -= this->PrecalculationShift;
  
      for (int i = firstComponent; i < LastComponent; ++i)
      {
        TmpNbrInteraction = this->NbrInteractionPerComponent[i];
        TmpIndexArray     = this->InteractionPerComponentIndex[i];
        TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
  
        Coefficient = vSource[k];
  
        for (j = 0; j < TmpNbrInteraction; ++j)
          vDestination[TmpIndexArray[j]] += TmpCoefficientArray[j] * Coefficient;
  
        vDestination[k++] += this->HamiltonianShift * Coefficient;
      }
    }
    else
    {
      if (this->DiskStorageFlag == false)
      {
        ParticleOnSphereWithNFlavor* TmpParticles = (ParticleOnSphereWithNFlavor*) this->Particles->Clone();
  
        int* TmpIndexArray;
        double* TmpCoefficientArray; 
        int j;
        int TmpNbrInteraction;
  
        firstComponent -= this->PrecalculationShift;
        LastComponent  -= this->PrecalculationShift;
  
        int Pos = firstComponent / this->FastMultiplicationStep; 
        int PosMod = firstComponent % this->FastMultiplicationStep;
  
        if (PosMod != 0)
        {
          ++Pos;
          PosMod = this->FastMultiplicationStep - PosMod;
        }
  
        int l = PosMod + firstComponent + this->PrecalculationShift;
  
        /* =============================
           Precalculated part
           ============================= */
  
        for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
        {
          TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
          TmpIndexArray = this->InteractionPerComponentIndex[Pos];
          TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
          Coefficient = vSource[l];
  
          for (j = 0; j < TmpNbrInteraction; ++j)
            vDestination[TmpIndexArray[j]] += TmpCoefficientArray[j] * Coefficient;
  
          vDestination[l] += this->HamiltonianShift * Coefficient;
  
          l += this->FastMultiplicationStep;
          ++Pos;
        }
  
        /* =============================
           Non-precalculated part
           ============================= */
  
        int Index;
        int* TmpIndices;
        double* TmpInteractionFactor;
        double Coefficient3;
  
        firstComponent += this->PrecalculationShift;
        LastComponent  += this->PrecalculationShift;
  
        for (l = 0; l < this->FastMultiplicationStep; ++l)
        if (PosMod != l)
        {
          for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
          {
            /* ===== INTRA ===== */
  
            for (int j = 0; j < this->NbrIntraSectorSums; ++j)
            {
              int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
              TmpIndices = this->IntraSectorIndicesPerSum[j];
  
              for (int i1 = 0; i1 < Lim; i1 += 2)
              {
                for (int sigma = 0; sigma < NbrFlavors; ++sigma)
                {
                  int channel = Channel(sigma, sigma);
  
                  Coefficient3 = TmpParticles->AsigmaAsigma( i, TmpIndices[i1], TmpIndices[i1+1], sigma, sigma);
  
                  if (Coefficient3 != 0.0)
                  {
                    TmpInteractionFactor = &(this->InteractionFactorsIntra[channel][j][(i1 * Lim) >> 2]);
                    Coefficient3 *= vSource[i];
  
                    for (int i2 = 0; i2 < Lim; i2 += 2)
                    {
                      Index = TmpParticles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2+1], sigma, sigma, Coefficient);
                      if (Index < Dim)
                        vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
  
                      ++TmpInteractionFactor;
                    }
                  }
                }
              }
            }
  
            /* ===== INTER ===== */
  
            for (int j = 0; j < this->NbrInterSectorSums; ++j)
            {
              int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
              TmpIndices = this->InterSectorIndicesPerSum[j];
  
              for (int i1 = 0; i1 < Lim; i1 += 2)
              {
                for (int s1 = 0; s1 < NbrFlavors; ++s1)
                for (int s2 = s1+1; s2 < NbrFlavors; ++s2)
                {
                  int channel = Channel(s1, s2);
  
                  Coefficient3 = TmpParticles->AsigmaAsigma( i, TmpIndices[i1], TmpIndices[i1+1], s1, s2);
  
                  if (Coefficient3 != 0.0)
                  {
                    TmpInteractionFactor = &(this->InteractionFactorsInter[channel][j][(i1 * Lim) >> 2]);
                    Coefficient3 *= vSource[i];
  
                    for (int i2 = 0; i2 < Lim; i2 += 2)
                    {
                      Index =
                        TmpParticles->AdsigmaAdsigma(
                          TmpIndices[i2],
                          TmpIndices[i2+1],
                          s1, s2,
                          Coefficient);
  
                      if (Index < Dim)
                        vDestination[Index] +=
                          Coefficient *
                          (*TmpInteractionFactor) *
                          Coefficient3;
  
                      ++TmpInteractionFactor;
                    }
                  }
                }
              }
            }    // Two-body inter interactions loop

            /* ===== ONE BODY (diag + tunneling) ===== */

            if (this->OneBodyInteractionFactors != 0)
            {
              double TmpDiagonal = 0.0;
              double Source = vSource[i];

              /* ---- diagonal ---- */
              for (int s = 0; s < NbrFlavors; ++s)
                for (int m = 0; m <= this->LzMax; ++m)
                  TmpDiagonal +=
                    this->OneBodyInteractionFactors[s][s][m] * TmpParticles->AdsigmaAsigma(i, m, s);

              vDestination[i] += (this->HamiltonianShift + TmpDiagonal) * Source;

              /* ---- tunneling ---- */
              int Index2;
              double Coeff2;

              for (int m = 0; m <= this->LzMax; ++m)
                for (int a = 0; a < NbrFlavors; ++a)
                  for (int b = 0; b < NbrFlavors; ++b)
                    if ((a != b) && (this->OneBodyInteractionFactors[a][b] != 0))
		      {
			Index2 =
			  TmpParticles->AdsigmaAsigma(i, m, m, a, b, Coeff2);

			if (Index2 < Dim)
			  vDestination[Index2] += Coeff2 * this->OneBodyInteractionFactors[a][b][m] * Source;
		      }
            }
            else
            {
              vDestination[i] += this->HamiltonianShift * vSource[i];
            }
          }
        }
  
        delete TmpParticles;
      }
      else
	    {
	      int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      int TmpNbrIteration = nbrComponent / this->BufferSize;
	      int* TmpIndexArray;
	      double* TmpCoefficientArray;
	      int TmpNbrInteraction;
	      int k = firstComponent;
	      int EffectiveHilbertSpaceDimension;
	      firstComponent -= this->PrecalculationShift;

	      ifstream File;
	      File.open(this->DiskStorageFileName, ios::binary | ios::in);
	      File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
	      long FileJump = 0;
	      for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
		FileJump += (long) this->NbrInteractionPerComponent[i];
	      FileJump *= sizeof(int);
	      long FileOffset = 0;
	      for (int i = this->DiskStorageStart; i < firstComponent; ++i)
		FileOffset += this->NbrInteractionPerComponent[i];
	      File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
	      FileJump += (sizeof(double) - sizeof(int)) * FileOffset;

	      for (int i = 0; i < TmpNbrIteration; ++i)
		{
		  int TmpPos = firstComponent;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < this->BufferSize; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);
		  
		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int l = 0; l < this->BufferSize; ++l)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      Coefficient = vSource[k];
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      vDestination[k] += this->HamiltonianShift * Coefficient;
		      ++k;
		      ++firstComponent;
		    }
		}

	      if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
		{
		  int TmpPos = firstComponent;
		  int Lim =  nbrComponent % this->BufferSize;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < Lim; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);

		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int i = 0; i < Lim; ++i)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      Coefficient = vSource[k];
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      vDestination[k] += this->HamiltonianShift * Coefficient;
		      ++k;
		      ++firstComponent;
		    }
		}

	      File.close();
	      delete[] BufferIndexArray;
	      delete[] BufferCoefficientArray;
	    }
    }
  }

  //cout << "Interaction : " << vDestination << endl;

  return vDestination;
}


RealVector* AbstractQHEOnSphereWithNFlavorHamiltonian::
LowLevelMultipleAddMultiply(RealVector* vSources,
                            RealVector* vDestinations,
                            int nbrVectors,
                            int firstComponent,
                            int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;

  if (this->FastMultiplicationFlag == false)
  {
    int Index;
    int* TmpIndices;
    double* TmpInteractionFactor;
    double Coefficient3;
    double* Coefficient2 = new double[nbrVectors];

    ParticleOnSphereWithNFlavor* TmpParticles = (ParticleOnSphereWithNFlavor*) this->Particles->Clone();

    for (int i = firstComponent; i < LastComponent; ++i)
    {
      // =========================
      // INTRA (σ = σ)
      // =========================
      for (int s = 0; s < NbrFlavors; ++s)
      {
        int channel = Channel(s,s);

        for (int j = 0; j < this->NbrIntraSectorSums; ++j)
        {
          int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
          TmpIndices = this->IntraSectorIndicesPerSum[j];

          for (int i1 = 0; i1 < Lim; i1 += 2)
          {
            Coefficient3 = TmpParticles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], s, s);

            if (Coefficient3 != 0.0)
            {
              TmpInteractionFactor = &(this->InteractionFactorsIntra[channel][j][(i1 * Lim) >> 2]);

              for (int p = 0; p < nbrVectors; ++p)
                Coefficient2[p] = Coefficient3 * vSources[p][i];

              for (int i2 = 0; i2 < Lim; i2 += 2)
              {
                Index =
                  TmpParticles->AdsigmaAdsigma(
                      TmpIndices[i2],
                      TmpIndices[i2 + 1], s, s, Coefficient);

                if (Index < Dim)
                  for (int p = 0; p < nbrVectors; ++p)
                    vDestinations[p][Index] +=
                        Coefficient *
                        (*TmpInteractionFactor) *
                        Coefficient2[p];

                ++TmpInteractionFactor;
              }
            }
          }
        }
      }

      // =========================
      // INTER (σ ≠ σ)
      // =========================
      for (int s1 = 0; s1 < NbrFlavors; ++s1)
      {
        for (int s2 = s1 + 1; s2 < NbrFlavors; ++s2)
        {
          int channel = Channel(s1,s2);

          for (int j = 0; j < this->NbrInterSectorSums; ++j)
          {
            int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
            TmpIndices = this->InterSectorIndicesPerSum[j];

            for (int i1 = 0; i1 < Lim; i1 += 2)
            {
              Coefficient3 = TmpParticles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], s1, s2);

              if (Coefficient3 != 0.0)
              {
                TmpInteractionFactor = &(this->InteractionFactorsInter[channel][j][(i1 * Lim) >> 2]);

                for (int p = 0; p < nbrVectors; ++p)
                  Coefficient2[p] = Coefficient3 * vSources[p][i];

                for (int i2 = 0; i2 < Lim; i2 += 2)
                {
                  Index =
                    TmpParticles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], s1, s2, Coefficient);

                  if (Index < Dim)
                    for (int p = 0; p < nbrVectors; ++p)
                      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];

                  ++TmpInteractionFactor;
                }
              }
            }
          }
        }
      }
    }

    // =========================
    // ONE BODY
    // =========================
    if (this->OneBodyInteractionFactors != 0)
    {
      int Dim = this->Particles->GetHilbertSpaceDimension();
      int Index;
      double Coefficient;
    
      for (int l = 0; l < nbrVectors; ++l)
      {
        RealVector& Src = vSources[l];
        RealVector& Dst = vDestinations[l];
    
        for (int i = firstComponent; i < LastComponent; ++i)
        {
          double TmpDiagonal = 0.0;
    
          for (int s = 0; s < NbrFlavors; ++s)
            for (int m = 0; m <= this->LzMax; ++m)
              TmpDiagonal +=
                this->OneBodyInteractionFactors[s][s][m] *
                TmpParticles->AdsigmaAsigma(i, m, s);
    
          Dst[i] += (this->HamiltonianShift + TmpDiagonal) * Src[i];
    
          for (int m = 0; m <= this->LzMax; ++m)
            for (int a = 0; a < NbrFlavors; ++a)
              for (int b = 0; b < NbrFlavors; ++b)
                if ((a != b) && (this->OneBodyInteractionFactors[a][b] != 0))
                {
                  Index = TmpParticles->AdsigmaAsigma(i, m, m, a, b, Coefficient);
    
                  if (Index < Dim)
                    Dst[Index] += Coefficient * this->OneBodyInteractionFactors[a][b][m] * Src[i];
                }
        }
      }
    }
    else
    {
      for (int l = 0; l < nbrVectors; ++l)
      {
        RealVector& Src = vSources[l];
        RealVector& Dst = vDestinations[l];

        for (int i = firstComponent; i < LastComponent; ++i)
          Dst[i] += this->HamiltonianShift * Src[i];
      }
    }

    delete[] Coefficient2;
    delete TmpParticles;
  }
  else
  {
    if (this->FastMultiplicationStep == 1)
    {
      int* TmpIndexArray;
      double* Coefficient2 = new double[nbrVectors];
      double* TmpCoefficientArray;
      int j;
      int Pos;
      int TmpNbrInteraction;
      int k = firstComponent;

      firstComponent -= this->PrecalculationShift;
      LastComponent -= this->PrecalculationShift;

      for (int i = firstComponent; i < LastComponent; ++i)
      {
        TmpNbrInteraction = this->NbrInteractionPerComponent[i];
        TmpIndexArray = this->InteractionPerComponentIndex[i];
        TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];

        for (int l = 0; l < nbrVectors; ++l)
        {
          Coefficient2[l] = vSources[l][k];
          vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
        }

        for (j = 0; j < TmpNbrInteraction; ++j)
        {
          Pos = TmpIndexArray[j];
          Coefficient = TmpCoefficientArray[j];

          for (int l = 0; l < nbrVectors; ++l)
            vDestinations[l][Pos] += Coefficient * Coefficient2[l];
        }

        ++k;
      }

      delete[] Coefficient2;
    }
    else
    {
      if (this->DiskStorageFlag == false)
        this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
      else
        this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  }

  return vDestinations;
}


RealVector* AbstractQHEOnSphereWithNFlavorHamiltonian::
LowLevelMultipleAddMultiplyPartialFastMultiply(
    RealVector* vSources, RealVector* vDestinations,
    int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;

  ParticleOnSphereWithNFlavor* TmpParticles =
      (ParticleOnSphereWithNFlavor*) this->Particles->Clone();

  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int j;
  int TmpNbrInteraction;

  firstComponent -= this->PrecalculationShift;
  LastComponent  -= this->PrecalculationShift;

  int Pos  = firstComponent / this->FastMultiplicationStep;
  int Pos2;
  double Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;

  double* Coefficient2 = new double[nbrVectors];
  int Dim = this->Particles->GetHilbertSpaceDimension();

  if (PosMod != 0)
  {
    ++Pos;
    PosMod = this->FastMultiplicationStep - PosMod;
  }

  int l = PosMod + firstComponent + this->PrecalculationShift;

  // ============================================================
  // PRECALCULATED PART
  // ============================================================

  for (int i = PosMod + firstComponent; i < LastComponent;
       i += this->FastMultiplicationStep)
  {
    TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
    TmpIndexArray     = this->InteractionPerComponentIndex[Pos];
    TmpCoefficientArray =
        this->InteractionPerComponentCoefficient[Pos];

    for (int k = 0; k < nbrVectors; ++k)
    {
      Coefficient2[k] = vSources[k][l];
      vDestinations[k][l] +=
          this->HamiltonianShift * Coefficient2[k];
    }

    for (j = 0; j < TmpNbrInteraction; ++j)
    {
      Pos2       = TmpIndexArray[j];
      Coefficient = TmpCoefficientArray[j];

      for (int k = 0; k < nbrVectors; ++k)
        vDestinations[k][Pos2] +=
            Coefficient * Coefficient2[k];
    }

    l += this->FastMultiplicationStep;
    ++Pos;
  }

  // ============================================================
  // NON-PRECALCULATED PART
  // ============================================================

  int Index;
  int* TmpIndices;
  double* TmpInteractionFactor;
  double Coefficient3;

  firstComponent += this->PrecalculationShift;
  LastComponent  += this->PrecalculationShift;

  for (l = 0; l < this->FastMultiplicationStep; ++l)
    if (PosMod != l)
    {
      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
      {
        // ============================
        // INTRA (σ = σ)
        // ============================

        for (int s = 0; s < NbrFlavors; ++s)
        {
          int channel = Channel(s,s);

          for (int j = 0; j < this->NbrIntraSectorSums; ++j)
          {
            int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
            TmpIndices = this->IntraSectorIndicesPerSum[j];

            for (int i1 = 0; i1 < Lim; i1 += 2)
            {
              Coefficient3 =
                TmpParticles->AsigmaAsigma(
                    i, TmpIndices[i1], TmpIndices[i1+1], s, s);

              if (Coefficient3 != 0.0)
              {
                TmpInteractionFactor =
                  &(this->InteractionFactorsIntra[channel][j]
                    [(i1 * Lim) >> 2]);

                for (int p = 0; p < nbrVectors; ++p)
                  Coefficient2[p] =
                      Coefficient3 * vSources[p][i];

                for (int i2 = 0; i2 < Lim; i2 += 2)
                {
                  Index =
                    TmpParticles->AdsigmaAdsigma(
                        TmpIndices[i2],
                        TmpIndices[i2+1],
                        s, s,
                        Coefficient);

                  if (Index < Dim)
                    for (int p = 0; p < nbrVectors; ++p)
                      vDestinations[p][Index] +=
                          Coefficient *
                          (*TmpInteractionFactor) *
                          Coefficient2[p];

                  ++TmpInteractionFactor;
                }
              }
            }
          }
        }

        // ============================
        // INTER (σ ≠ σ)
        // ============================

        for (int s1 = 0; s1 < NbrFlavors; ++s1)
        {
          for (int s2 = s1+1; s2 < NbrFlavors; ++s2)
          {
            int channel = Channel(s1,s2);

            for (int j = 0; j < this->NbrInterSectorSums; ++j)
            {
              int Lim = 2 *
                this->NbrInterSectorIndicesPerSum[j];
              TmpIndices =
                this->InterSectorIndicesPerSum[j];

              for (int i1 = 0; i1 < Lim; i1 += 2)
              {
                Coefficient3 =
                  TmpParticles->AsigmaAsigma(
                      i,
                      TmpIndices[i1],
                      TmpIndices[i1+1],
                      s1, s2);

                if (Coefficient3 != 0.0)
                {
                  TmpInteractionFactor =
                    &(this->InteractionFactorsInter[channel][j]
                      [(i1 * Lim) >> 2]);

                  for (int p = 0; p < nbrVectors; ++p)
                    Coefficient2[p] =
                        Coefficient3 * vSources[p][i];

                  for (int i2 = 0; i2 < Lim; i2 += 2)
                  {
                    Index =
                      TmpParticles->AdsigmaAdsigma(
                          TmpIndices[i2],
                          TmpIndices[i2+1],
                          s1, s2,
                          Coefficient);

                    if (Index < Dim)
                      for (int p = 0; p < nbrVectors; ++p)
                        vDestinations[p][Index] +=
                            Coefficient *
                            (*TmpInteractionFactor) *
                            Coefficient2[p];

                    ++TmpInteractionFactor;
                  }
                }
              }
            }
          }
        } // Two-body inter interactions loop


        // ============================ 
        // ONE BODY
        // ============================ 
        if (this->OneBodyInteractionFactors != 0)
        {
          double TmpDiagonal = 0.0;
          for (int s = 0; s < NbrFlavors; ++s)
            for (int m = 0; m <= this->LzMax; ++m)
              TmpDiagonal +=
                this->OneBodyInteractionFactors[s][s][m] *
                TmpParticles->AdsigmaAsigma(i, m, s);

          for (int p = 0; p < nbrVectors; ++p)
            vDestinations[p][i] +=
                (this->HamiltonianShift + TmpDiagonal) * vSources[p][i];

          int Index;
          double Coefficient;
          for (int m = 0; m <= this->LzMax; ++m)
            for (int a = 0; a < NbrFlavors; ++a)
              for (int b = 0; b < NbrFlavors; ++b)
                if ((a != b) && (this->OneBodyInteractionFactors[a][b] != 0))
		  {
		    Index =
		      TmpParticles->AdsigmaAsigma(
						  i, m, m, a, b, Coefficient);
		    
		    if (Index < Dim)
		      for (int p = 0; p < nbrVectors; ++p)
			vDestinations[p][Index] +=
                          Coefficient * this->OneBodyInteractionFactors[a][b][m] * vSources[p][i];
		  }
        }   // One-body loop
      }    // Non-precalculated part loop
    }   // Fast multiplication step loop

  delete[] Coefficient2;
  delete TmpParticles;

  return vDestinations;
}


long AbstractQHEOnSphereWithNFlavorHamiltonian::
PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  long Memory = 0;

  int* TmpIndices;

  ParticleOnSphereWithNFlavor* TmpParticles =
      (ParticleOnSphereWithNFlavor*) this->Particles->Clone();

  int LastComponent = lastComponent + firstComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();

  for (int i = firstComponent; i < LastComponent; ++i)
  {
    // ============================
    // INTRA (σ = σ)
    // ============================
    for (int s = 0; s < NbrFlavors; ++s)
    {
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
      {
        int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
        TmpIndices = this->IntraSectorIndicesPerSum[j];

        for (int i1 = 0; i1 < Lim; i1 += 2)
        {
          Coefficient2 =
            TmpParticles->AsigmaAsigma(
                i,
                TmpIndices[i1],
                TmpIndices[i1 + 1],
                s, s);

          if (Coefficient2 != 0.0)
          {
            for (int i2 = 0; i2 < Lim; i2 += 2)
            {
              Index =
                TmpParticles->AdsigmaAdsigma(
                    TmpIndices[i2],
                    TmpIndices[i2 + 1],
                    s, s,
                    Coefficient);

              if (Index < Dim)
              {
                ++Memory;
                ++this->NbrInteractionPerComponent[
                    i - this->PrecalculationShift];
              }
            }
          }
        }
      }
    }

    // ============================
    // INTER (σ ≠ σ)
    // ============================
    for (int s1 = 0; s1 < NbrFlavors; ++s1)
    {
      for (int s2 = s1 + 1; s2 < NbrFlavors; ++s2)
      {
        for (int j = 0; j < this->NbrInterSectorSums; ++j)
        {
          int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
          TmpIndices = this->InterSectorIndicesPerSum[j];

          for (int i1 = 0; i1 < Lim; i1 += 2)
          {
            Coefficient2 =
              TmpParticles->AsigmaAsigma(
                  i,
                  TmpIndices[i1],
                  TmpIndices[i1 + 1],
                  s1, s2);

            if (Coefficient2 != 0.0)
            {
              for (int i2 = 0; i2 < Lim; i2 += 2)
              {
                Index =
                  TmpParticles->AdsigmaAdsigma(
                      TmpIndices[i2],
                      TmpIndices[i2 + 1],
                      s1, s2,
                      Coefficient);

                if (Index < Dim)
                {
                  ++Memory;
                  ++this->NbrInteractionPerComponent[
                      i - this->PrecalculationShift];
                }
              }
            }
          }
        }
      }
    }

    // ============================
    // ONE BODY
    // ============================
    if (this->OneBodyInteractionFactors != 0)
    {
      /* diagonal */
      ++Memory;
      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
    
      /* tunneling */
      for (int m = 0; m <= this->LzMax; ++m)
        for (int a = 0; a < NbrFlavors; ++a)
          for (int b = 0; b < NbrFlavors; ++b)
            if ((a != b) && (this->OneBodyInteractionFactors[a][b] != 0))
	      {
		Index = TmpParticles->AdsigmaAsigma(i, m, m, a, b, Coefficient);
		
		if (Index < Dim)
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
	      }
    }
  
  
  }

  delete TmpParticles;
  return Memory;
}


void AbstractQHEOnSphereWithNFlavorHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;

  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);

  int EffectiveHilbertSpaceDimension =
      ((int)(MaxIndex - MinIndex)) + 1;

  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;

  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;

  gettimeofday(&(TotalStartingTime2), 0);
  cout << "start" << endl;

  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;

  if ((ReducedSpaceDimension *
       this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;

  int* TmpNbrInteractionPerComponent = this->NbrInteractionPerComponent;
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;


  
  //this->PartialFastMultiplicationMemory( this->PrecalculationShift, EffectiveHilbertSpaceDimension);
  
  this->InteractionPerComponentIndex = new int*[ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double*[ReducedSpaceDimension];

  int Dim = this->Particles->GetHilbertSpaceDimension();
  double* TmpInteractionFactor;
  int* TmpIndices;

  ParticleOnSphereWithNFlavor* TmpParticles =
      (ParticleOnSphereWithNFlavor*) this->Particles->Clone();

  int TotalPos = 0;

  for (int i = 0;
       i < EffectiveHilbertSpaceDimension;
       i += this->FastMultiplicationStep)
  {
    this->NbrInteractionPerComponent[TotalPos] = TmpNbrInteractionPerComponent[i];
    this->InteractionPerComponentIndex[TotalPos] = new int[this->NbrInteractionPerComponent[TotalPos]];

    this->InteractionPerComponentCoefficient[TotalPos] = new double[this->NbrInteractionPerComponent[TotalPos]];

    // cout << "i = " << i << " " << this->PrecalculationShift << endl;
    // cout << "TotalPos = " << TotalPos << " " << this->NbrInteractionPerComponent[TotalPos] << endl;

    TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
    TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];

    Pos = 0;

    // ============================
    // INTRA (σ = σ)
    // ============================

    for (int s = 0; s < NbrFlavors; ++s)
    {
      int channel = Channel(s,s);

      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
      {
        int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
        TmpIndices = this->IntraSectorIndicesPerSum[j];

        for (int i1 = 0; i1 < Lim; i1 += 2)
        {
          Coefficient2 =
            TmpParticles->AsigmaAsigma(
                i + this->PrecalculationShift,
                TmpIndices[i1],
                TmpIndices[i1+1],
                s, s);

          if (Coefficient2 != 0.0)
          {
            TmpInteractionFactor =
              &(this->InteractionFactorsIntra[channel][j]
                [(i1 * Lim) >> 2]);

            for (int i2 = 0; i2 < Lim; i2 += 2)
            {
              Index =
                TmpParticles->AdsigmaAdsigma(
                    TmpIndices[i2],
                    TmpIndices[i2+1],
                    s, s,
                    Coefficient);

              if (Index < Dim)
              {
                TmpIndexArray[Pos] = Index;
                TmpCoefficientArray[Pos] =
                    Coefficient *
                    Coefficient2 *
                    (*TmpInteractionFactor);
                ++Pos;
              }
              ++TmpInteractionFactor;
            }
          }
        }
      }
    }

    // ============================
    // INTER (σ ≠ σ)
    // ============================

    for (int s1 = 0; s1 < NbrFlavors; ++s1)
    {
      for (int s2 = s1+1; s2 < NbrFlavors; ++s2)
      {
        int channel = Channel(s1,s2);

        for (int j = 0; j < this->NbrInterSectorSums; ++j)
        {
          int Lim =
            2 * this->NbrInterSectorIndicesPerSum[j];

          TmpIndices =
            this->InterSectorIndicesPerSum[j];

          for (int i1 = 0; i1 < Lim; i1 += 2)
          {
            Coefficient2 =
              TmpParticles->AsigmaAsigma(
                  i + this->PrecalculationShift,
                  TmpIndices[i1],
                  TmpIndices[i1+1],
                  s1, s2);

            if (Coefficient2 != 0.0)
            {
              TmpInteractionFactor =
                &(this->InteractionFactorsInter[channel][j]
                  [(i1 * Lim) >> 2]);

              for (int i2 = 0; i2 < Lim; i2 += 2)
              {
                Index =
                  TmpParticles->AdsigmaAdsigma(
                      TmpIndices[i2],
                      TmpIndices[i2+1],
                      s1, s2,
                      Coefficient);

                if (Index < Dim)
                {
                  TmpIndexArray[Pos] = Index;
                  TmpCoefficientArray[Pos] =
                      Coefficient *
                      Coefficient2 *
                      (*TmpInteractionFactor);
                  ++Pos;
                }
                ++TmpInteractionFactor;
              }
            }
          }
        }
      }
    }

    // ============================
    // ONE BODY (diag + tunneling)
    // ============================

    if (this->OneBodyInteractionFactors != 0)
    {
      double TmpDiagonal = 0.0;

      for (int s = 0; s < NbrFlavors; ++s)
	{
	  if (this->OneBodyInteractionFactors[s][s] !=0)
	    {
	      for (int m = 0; m <= this->LzMax; ++m)
		TmpDiagonal += this->OneBodyInteractionFactors[s][s][m] *
		  TmpParticles->AdsigmaAsigma(i + this->PrecalculationShift, m, s);
	    }
	}
      TmpIndexArray[Pos] = i + this->PrecalculationShift;
      TmpCoefficientArray[Pos] = TmpDiagonal;   // + this->HamiltonianShift
      ++Pos;

      int Index;
      double Coefficient;

      for (int m = 0; m <= this->LzMax; ++m)
        for (int a = 0; a < NbrFlavors; ++a)
          for (int b = 0; b < NbrFlavors; ++b)
            if ((a != b) && (this->OneBodyInteractionFactors[a][b] != 0))
            {
              Index =
                TmpParticles->AdsigmaAsigma(i + this->PrecalculationShift, m, m, a, b, Coefficient);

              if (Index < Dim)
              {
                TmpIndexArray[Pos] = Index;
                TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[a][b][m];    // + this->HamiltonianShift
                ++Pos;
              }
            }
    }

    ++TotalPos;
  }
  delete[] TmpNbrInteractionPerComponent;
  
  this->FastMultiplicationFlag = true;
  gettimeofday(&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;

  Dt2 = (double)(TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}



void AbstractQHEOnSphereWithNFlavorHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
}


void AbstractQHEOnSphereWithNFlavorHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  if (this->FastMultiplicationStep == 1)
    {
      this->DiskStorageFlag = false;
      this->DiskStorageFileName = 0;
      this->EnableFastMultiplication();
      return;
    }

  this->DiskStorageFlag = true;
  this->DiskStorageFileName = new char [strlen(fileName) + 8];
  sprintf(this->DiskStorageFileName, "%s.ham", fileName);

  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int)(MaxIndex - MinIndex)) + 1;

  this->DiskStorageStart = (int)MinIndex;
  int DiskStorageEnd = 1 + (int)MaxIndex;

  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;

  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;

  gettimeofday(&(TotalStartingTime2), 0);
  cout << "start" << endl;

  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  this->MaxNbrInteractionPerComponent = 0;

  int TotalPos = 0;

  ofstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::out);

  File.write((char*) &(EffectiveHilbertSpaceDimension), sizeof(int));
  File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
  File.write((char*) this->NbrInteractionPerComponent,
             sizeof(int) * EffectiveHilbertSpaceDimension);

  long FileJump = 0;

  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      FileJump += (long)this->NbrInteractionPerComponent[i];

      if (this->MaxNbrInteractionPerComponent <
          this->NbrInteractionPerComponent[i])
        this->MaxNbrInteractionPerComponent =
            this->NbrInteractionPerComponent[i];
    }

  FileJump *= sizeof(int);

  TmpIndexArray = new int[this->MaxNbrInteractionPerComponent];
  TmpCoefficientArray = new double[this->MaxNbrInteractionPerComponent];

  for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
    {
      if (this->NbrInteractionPerComponent[TotalPos] > 0)
        {
          Pos = 0;

          File.write((char*) TmpIndexArray,
                     sizeof(int) *
                     this->NbrInteractionPerComponent[TotalPos]);

          FileJump -= sizeof(int) *
                      this->NbrInteractionPerComponent[TotalPos];

          File.seekp(FileJump, ios::cur);

          File.write((char*) TmpCoefficientArray,
                     sizeof(double) *
                     this->NbrInteractionPerComponent[TotalPos]);

          FileJump += sizeof(double) *
                      this->NbrInteractionPerComponent[TotalPos];

          File.seekp(-FileJump, ios::cur);
        }

      ++TotalPos;
    }

  delete[] TmpIndexArray;
  delete[] TmpCoefficientArray;
  File.close();

  this->FastMultiplicationFlag = true;

  this->BufferSize =
    this->Memory /
    ((this->MaxNbrInteractionPerComponent *
     (sizeof(int) + sizeof(double)))
     + sizeof(int*) + sizeof(double*));

  gettimeofday(&(TotalEndingTime2), 0);

  cout << "------------------------------------------------------------------"
       << endl << endl;

  Dt2 =
    (double)(TotalEndingTime2.tv_sec -
             TotalStartingTime2.tv_sec)
    +
    ((TotalEndingTime2.tv_usec -
      TotalStartingTime2.tv_usec) / 1000000.0);

  cout << "time = " << Dt2 << endl;
}


