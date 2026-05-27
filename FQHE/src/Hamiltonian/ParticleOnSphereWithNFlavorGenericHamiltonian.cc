////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of N-flavor Hamiltonian for particles on a sphere            //
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

#include "Hamiltonian/ParticleOnSphereWithNFlavorGenericHamiltonian.h"
#include "HilbertSpace/FermionOnSphereWithNFlavor.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

using std::cout;
using std::endl;


ParticleOnSphereWithNFlavorGenericHamiltonian::ParticleOnSphereWithNFlavorGenericHamiltonian(ParticleOnSphereWithNFlavor* particles, int nbrParticles, int lzmax, double** pseudoPotential, double*** oneBodyPotential, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->NbrFlavors = particles->GetNbrFlavors();

  int NChannels = this->NbrFlavors * (this->NbrFlavors + 1) / 2;

  this->PseudoPotentials = new double* [NChannels];
  for (int j = 0; j < NChannels; ++j)
  {
    this->PseudoPotentials[j] = new double [this->NbrLzValue];
    for (int i = 0; i < this->NbrLzValue; ++i)
    this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
  }

  this->OneBodyInteractionFactors = 0;

  if (oneBodyPotential)
  {
    this->OneBodyTermFlag = true;

    this->OneBodyInteractionFactors = new double** [this->NbrFlavors];
    for (int a = 0; a < this->NbrFlavors; ++a)
    {
      this->OneBodyInteractionFactors[a] = new double* [this->NbrFlavors];

      for (int b = 0; b < this->NbrFlavors; ++b)
      {
        this->OneBodyInteractionFactors[a][b] = new double[this->NbrLzValue];

        for (int m = 0; m < this->NbrLzValue; ++m)
          this->OneBodyInteractionFactors[a][b][m] = oneBodyPotential[a][b][m];
      }
    }
  }

  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;

  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;

  if (precalculationFileName == 0)
  {
    if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);

      if (TmpMemory < 1024)
      cout  << "fast = " <<  TmpMemory << "b ";
      else
      if (TmpMemory < (1 << 20))
      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
      else
      if (TmpMemory < (1 << 30))
      cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
      else
      {
        cout  << "fast = " << (TmpMemory >> 30) << ".";
        TmpMemory -= ((TmpMemory >> 30) << 30);
        TmpMemory *= 100l;
        TmpMemory >>= 30;
        if (TmpMemory < 10l)
        cout << "0";
        cout  << TmpMemory << " Gb ";
      }

      if (this->DiskStorageFlag == false)
      this->EnableFastMultiplication();
      else
      {
        char* TmpFileName = this->Architecture->GetTemporaryFileName();
        this->EnableFastMultiplicationWithDiskStorage(TmpFileName);
        delete[] TmpFileName;
      }
    }
  }
  else
  this->LoadPrecalculation(precalculationFileName);
}


ParticleOnSphereWithNFlavorGenericHamiltonian::~ParticleOnSphereWithNFlavorGenericHamiltonian() 
{
  int NChannels = this->NbrFlavors * (this->NbrFlavors + 1) / 2;

  for (int j = 0; j < NChannels; ++j)
    delete[] this->PseudoPotentials[j];

  delete[] this->PseudoPotentials;
}


void ParticleOnSphereWithNFlavorGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnSphereWithNFlavor*) hilbertSpace;
  this->EvaluateInteractionFactors();
}


AbstractHilbertSpace* ParticleOnSphereWithNFlavorGenericHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}


int ParticleOnSphereWithNFlavorGenericHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}


void ParticleOnSphereWithNFlavorGenericHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}


Complex ParticleOnSphereWithNFlavorGenericHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  cout << "Testing matrix element H(0,0) = " << this->MatrixElement(V1, V2) << endl;
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}

Complex ParticleOnSphereWithNFlavorGenericHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  cout << "Testing matrix element H(0,0) = " << this->MatrixElement(V1, V2) << endl;
  return Complex();
}


List<Matrix*> ParticleOnSphereWithNFlavorGenericHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}


List<Matrix*> ParticleOnSphereWithNFlavorGenericHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}



void ParticleOnSphereWithNFlavorGenericHamiltonian::EvaluateInteractionFactors()
{
  //cout << "Entering EvaluateInteractionFactors()" << endl;
  int Lim;
  int Min;
  int Pos = 0;
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  int J = 2 * this->LzMax - 2;
  int m4;
  double ClebschCoef;
  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;

  this->NbrInterSectorSums = 2 * this->LzMax + 1;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];

  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;

  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrInterSectorIndicesPerSum[m1 + m2];

  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];

  for (int i = 0; i < this->NbrInterSectorSums; ++i)
  {
    this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];
    this->NbrInterSectorIndicesPerSum[i] = 0;
  }

  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
    {
      this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
      this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
      ++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
    }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
  {
    this->NbrIntraSectorSums = 2 * this->LzMax - 1;
    this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];

    for (int i = 0; i < this->NbrIntraSectorSums; ++i)
      this->NbrIntraSectorIndicesPerSum[i] = 0;

    for (int m1 = 0; m1 < this->LzMax; ++m1)
      for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
        ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];

      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];

    for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    {
      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];
      this->NbrIntraSectorIndicesPerSum[i] = 0;
    }

    for (int m1 = 0; m1 < this->LzMax; ++m1)
      for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
      {
        this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
        this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
        ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
      }

    int N = this->NbrFlavors;
    int NChannels = N * (N + 1) / 2;

    this->InteractionFactorsIntra = new double** [NChannels];
    this->InteractionFactorsInter = new double** [NChannels];

    for (int c = 0; c < NChannels; ++c) this->InteractionFactorsIntra[c] = 0;
    for (int c = 0; c < NChannels; ++c) this->InteractionFactorsInter[c] = 0;

    for (int a = 0; a < N; ++a)
      for (int b = a; b < N; ++b)
      {
        int channel = Channel(a,b);

        if (a == b)
        {
          this->InteractionFactorsIntra[channel] = new double* [this->NbrIntraSectorSums];

          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
          {
            this->InteractionFactorsIntra[channel][i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];

            int Index = 0;

            for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
            {
              int m1 = (this->IntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
              int m2 = (this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;

              for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
              {
                int m3 = (this->IntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
                int m4 = (this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;

                Clebsch.InitializeCoefficientIterator(m1, m2);

                this->InteractionFactorsIntra[channel][i][Index] = 0.0;

                while (Clebsch.Iterate(J, ClebschCoef))
                  if (((J >> 1) & 1) == Sign)
                  {
                    TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);

                    this->InteractionFactorsIntra[channel][i][Index] += this->PseudoPotentials[channel][J >> 1] * TmpCoefficient;
                  }

                this->InteractionFactorsIntra[channel][i][Index] *= -4.0;
                //if (a==0 && b==0)
                //cout << a<<b<< ": Channel = " << Channel(a,b) << " = " << channel << ":: (j1, j2, Index, TotalNbrInteractionFactors) = " << j1<<j2<< Index << TotalNbrInteractionFactors << "] = " << this->InteractionFactorsIntra[channel][i][Index] << endl;

                ++TotalNbrInteractionFactors;
                ++Index;
              }
            }
          }
        }
        else
        {
          this->InteractionFactorsInter[channel] = new double* [this->NbrInterSectorSums];

          for (int i = 0; i < this->NbrInterSectorSums; ++i)
          {
            this->InteractionFactorsInter[channel][i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];

            int Index = 0;

            for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
            {
              int m1 = (this->InterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
              int m2 = (this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;

              for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
              {
                int m3 = (this->InterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
                int m4 = (this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;

                Clebsch.InitializeCoefficientIterator(m1, m2);

                this->InteractionFactorsInter[channel][i][Index] = 0.0;

                while (Clebsch.Iterate(J, ClebschCoef))
                {
                  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
                  this->InteractionFactorsInter[channel][i][Index] += this->PseudoPotentials[channel][J >> 1] * TmpCoefficient;
                }

                this->InteractionFactorsInter[channel][i][Index] *= -2.0;
                //if (a==0 && b==1) cout << a<<b<< ": Channel = " << Channel(a,b) << " = " << channel << ":: (j1, j2, Index, TotalNbrInteractionFactors) = " << j1<<j2<< Index << TotalNbrInteractionFactors << "] = " << this->InteractionFactorsInter[channel][i][Index] << endl;

                ++TotalNbrInteractionFactors;
                ++Index;
              }
            }
          }
        }
      
      }
  }

  else
  {
  }


  cout << "nbr interaction = "
       << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}




