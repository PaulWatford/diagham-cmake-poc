////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//        class of a two body interaction projected onto a single band        //
//        from an ASCII file providing the three body matrix elements         //
//                                                                            //
//                        last modification : 07/05/2024                      //
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
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// matrixElementsThreeBodyInteractionFile = name of the ASCII file containing the matrix element for the generic three body interaction term
// threeBodyInteractionRescalingFactor = global rescaling factor for the three-body interaction term
// matrixElementsTwoBodyInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
// twoBodyInteractionRescalingFactor = global rescaling factor for the two-body interaction term
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// spinFlag = include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,	
																 char* matrixElementsThreeBodyInteractionFile,
																 double threeBodyInteractionRescalingFactor,
																 char* matrixElementsTwoBodyInteractionFile,
																 double twoBodyInteractionRescalingFactor, 
																 Abstract2DTightBindingModel* tightBindingModel ,bool flatBandFlag,
																 bool spinFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->ThreeBodyInteractionRescalingFactor = threeBodyInteractionRescalingFactor;
  this->TwoBodyInteractionRescalingFactor = twoBodyInteractionRescalingFactor;
  this->AdditionalSpinFlag = spinFlag;
  if (this->AdditionalSpinFlag == true)
    {
      this->NbrInternalIndices = 2;
    }
  else
    {
      this->NbrInternalIndices = 1;
    }
  this->MatrixElementsThreeBodyInteractionFile = new char[strlen(matrixElementsThreeBodyInteractionFile) + 1];
  strcpy(this->MatrixElementsThreeBodyInteractionFile, matrixElementsThreeBodyInteractionFile);
  if (matrixElementsTwoBodyInteractionFile == 0)
    {
      this->MatrixElementsTwoBodyInteractionFile = 0;
    }
  else
    {
      this->MatrixElementsTwoBodyInteractionFile = new char[strlen(matrixElementsTwoBodyInteractionFile) + 1];
      strcpy(this->MatrixElementsTwoBodyInteractionFile, matrixElementsTwoBodyInteractionFile);
    }
  
  this->InitializeNBodyInteraction();
  
  this->Architecture = architecture;
  this->Memory = memory;
  // this->OneBodyInteractionFactorsupup = 0;
  // this->OneBodyInteractionFactorsdowndown = 0;
  // this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;//false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
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
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::~ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian()
{
}

// initialize the n-body interaction terms
//

void ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::InitializeNBodyInteraction()
{
  this->MaxNBody = 3;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NbrSpinSectors = new int [this->MaxNBody + 1];
  this->NBodySign = new double*[this->MaxNBody + 1];
  this->SpinIndices = new int** [this->MaxNBody + 1];
  this->SpinIndicesShort = new int* [this->MaxNBody + 1];
  this->NbrNBodySpinMomentumSectorSum = new int* [this->MaxNBody + 1];
  this->NbrNBodySpinMomentumSectorIndicesPerSum = new int**[this->MaxNBody + 1];
  this->NBodySpinMomentumSectorIndicesPerSum = new int***[this->MaxNBody + 1];
  this->NBodyInteractionFactors = new Complex***[this->MaxNBody + 1];
  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->NBodyFlags[k] = false;
      this->NbrSpinSectors[k] = 0;
      this->NBodySign[k] = 0;
      this->SpinIndices[k] = 0;
      this->SpinIndicesShort[k] = 0;
      this->NBodyInteractionFactors[k] = 0;
    }

  if (this->AdditionalSpinFlag == false)
    {
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  this->NBodySign[3] = new double[1];
	  this->NBodySign[3][0] = -1.0;
	}
      this->NbrSpinSectors[3] = 1;
      this->SpinIndices[3] = new int*[this->NbrSpinSectors[3]];
      this->SpinIndices[3][0] = new int[6];
      this->SpinIndices[3][0][0] = 0;
      this->SpinIndices[3][0][1] = 0;
      this->SpinIndices[3][0][2] = 0;
      this->SpinIndices[3][0][3] = 0;
      this->SpinIndices[3][0][4] = 0;
      this->SpinIndices[3][0][5] = 0;
      this->SpinIndicesShort[3] = new int[this->NbrSpinSectors[3]];
      this->SpinIndicesShort[3][0] = 0x0;
    }
  else
    {
       if (this->AdditionalSpinFlag == false)
	{
	   this->NBodySign[3] = new double[4];
	   this->NBodySign[3][0] = -1.0;
	   this->NBodySign[3][1] = -1.0;
	   this->NBodySign[3][2] = 1.0;
	   this->NBodySign[3][3] = 1.0;
	 }
      this->NbrSpinSectors[3] = 4;
      this->SpinIndices[3] = new int*[this->NbrSpinSectors[3]];
      this->SpinIndices[3][0] = new int[6];
      this->SpinIndices[3][1] = new int[6];
      this->SpinIndices[3][2] = new int[6];
      this->SpinIndices[3][3] = new int[6];

      this->SpinIndices[3][0][0] = 0;
      this->SpinIndices[3][0][1] = 0;
      this->SpinIndices[3][0][2] = 0;
      this->SpinIndices[3][0][3] = 0;
      this->SpinIndices[3][0][4] = 0;
      this->SpinIndices[3][0][5] = 0;
      this->SpinIndices[3][1][0] = 1;
      this->SpinIndices[3][1][1] = 1;
      this->SpinIndices[3][1][2] = 1;
      this->SpinIndices[3][1][3] = 1;
      this->SpinIndices[3][1][4] = 1;
      this->SpinIndices[3][1][5] = 1;
      
      this->SpinIndices[3][2][0] = 0;
      this->SpinIndices[3][2][1] = 1;
      this->SpinIndices[3][2][2] = 1;
      this->SpinIndices[3][2][3] = 0;
      this->SpinIndices[3][2][4] = 1;
      this->SpinIndices[3][2][5] = 1;
      
      this->SpinIndices[3][3][0] = 1;
      this->SpinIndices[3][3][1] = 0;
      this->SpinIndices[3][3][2] = 0;
      this->SpinIndices[3][3][3] = 1;
      this->SpinIndices[3][3][4] = 0;
      this->SpinIndices[3][3][5] = 0;
      this->SpinIndicesShort[3] = new int[this->NbrSpinSectors[3]];
      this->SpinIndicesShort[3][0] = 0x0;
      this->SpinIndicesShort[3][1] = 0x7 | (0x7<<3);
      this->SpinIndicesShort[3][2] = 0x4 | (0x4<<3);
      this->SpinIndicesShort[3][3] = 0x3 | (0x3<<3);
    }
  
  this->NBodyFlags[3] = true;
  this->NbrNBodySpinMomentumSectorSum[3] = new int[this->NbrSpinSectors[3]];
  this->NbrNBodySpinMomentumSectorIndicesPerSum[3] = new int* [this->NbrSpinSectors[3]];
  this->NBodySpinMomentumSectorIndicesPerSum[3] = new int** [this->NbrSpinSectors[3]];
  this->NBodyInteractionFactors[3] = new Complex**[this->NbrSpinSectors[3]];
}

// evaluate all interaction factors
//   

void ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::EvaluateInteractionFactors()
{
  this->EvaluateThreeBodyInteractionFactors();
  this->EvaluateTwoBodyInteractionFactors();
  this->EvaluateOneBodyInteractionFactors(); 
}

// evaluate the three-body interaction factors 
//   

void ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::EvaluateThreeBodyInteractionFactors()
{
  for (int s = 0; s < this->NbrSpinSectors[3]; ++s)
    {
      this->NbrNBodySpinMomentumSectorSum[3][s] = this->NbrSiteX * this->NbrSiteY;
      this->NbrNBodySpinMomentumSectorIndicesPerSum[3][s] = new int[this->NbrNBodySpinMomentumSectorSum[3][s]];
      this->NBodySpinMomentumSectorIndicesPerSum[3][s] = new int*[this->NbrNBodySpinMomentumSectorSum[3][s]];
    }


  Complex* TmpArrayMatrixElements = 0;
  int* TmpArraySigma1 = 0;
  int* TmpArraySigma2 = 0;
  int* TmpArraySigma3 = 0;
  int* TmpArraySigma4 = 0;
  int* TmpArraySigma5 = 0;
  int* TmpArraySigma6 = 0;
  int* TmpArrayKx1 = 0;
  int* TmpArrayKy1 = 0;
  int* TmpArrayKx2 = 0;
  int* TmpArrayKy2 = 0;
  int* TmpArrayKx3 = 0;
  int* TmpArrayKy3 = 0;
  int* TmpArrayKx4 = 0;
  int* TmpArrayKy4 = 0;
  int* TmpArrayKx5 = 0;
  int* TmpArrayKy5 = 0;
  int* TmpArrayKx6 = 0;
  int* TmpArrayKy6 = 0;
  int TmpNbrThreeBodyMatrixElements = this->ProcessThreeBodyMatrixElements(TmpArraySigma1, TmpArraySigma2, TmpArraySigma3, TmpArraySigma4, TmpArraySigma5, TmpArraySigma6,
									   TmpArrayKx1, TmpArrayKy1, TmpArrayKx2, TmpArrayKy2, TmpArrayKx3, TmpArrayKy3,
									   TmpArrayKx4, TmpArrayKy4, TmpArrayKx5, TmpArrayKy5, TmpArrayKx6, TmpArrayKy6, TmpArrayMatrixElements);
  int* TmpLinearizedSumK = new int[TmpNbrThreeBodyMatrixElements];
  int* TmpLinearizedK1 = new int[TmpNbrThreeBodyMatrixElements];
  int* TmpLinearizedK2 = new int[TmpNbrThreeBodyMatrixElements];
  int* TmpLinearizedK3 = new int[TmpNbrThreeBodyMatrixElements];
  int* TmpLinearizedK4 = new int[TmpNbrThreeBodyMatrixElements];
  int* TmpLinearizedK5 = new int[TmpNbrThreeBodyMatrixElements];
  int* TmpLinearizedK6 = new int[TmpNbrThreeBodyMatrixElements];
  for (int i = 0; i < TmpNbrThreeBodyMatrixElements; ++i)
    {
      TmpLinearizedSumK[i] = this->TightBindingModel->GetLinearizedMomentumIndex((TmpArrayKx1[i] + TmpArrayKx2[i] + TmpArrayKx3[i]) % this->NbrSiteX,
										 (TmpArrayKy1[i] + TmpArrayKy2[i] + TmpArrayKy3[i]) % this->NbrSiteY);
      TmpLinearizedK1[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpArrayKx1[i], TmpArrayKy1[i]);
      TmpLinearizedK2[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpArrayKx2[i], TmpArrayKy2[i]);
      TmpLinearizedK3[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpArrayKx3[i], TmpArrayKy3[i]);
      TmpLinearizedK4[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpArrayKx4[i], TmpArrayKy4[i]);
      TmpLinearizedK5[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpArrayKx5[i], TmpArrayKy5[i]);
      TmpLinearizedK6[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpArrayKx6[i], TmpArrayKy6[i]);
      TmpArrayMatrixElements[i] *= this->ThreeBodyInteractionRescalingFactor;
    }

  
  int** Permutations = 0; 
  double* PermutationSign = 0; 
  int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign, 3);
  
  if (this->AdditionalSpinFlag == false)
    {
    }
  else
    {
      cout << "spinful case not yet implemented" << endl;
      exit(0);
    }
}

// evaluate  the two-body interaction factors
//   

void ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::EvaluateTwoBodyInteractionFactors()
{
  long TotalNbrInteractionFactors = 0l;

  int* TmpSigma1 = 0;
  int* TmpSigma2 = 0;
  int* TmpSigma3 = 0;
  int* TmpSigma4 = 0;
  int* TmpKx1 = 0;
  int* TmpKx2 = 0;
  int* TmpKx3 = 0;
  int* TmpKx4 = 0;
  int* TmpKy1 = 0;
  int* TmpKy2 = 0;
  int* TmpKy3 = 0;
  int* TmpKy4 = 0;
  Complex* TmpMatrixElements = 0;

  int TmpNbrTwoBodyMatrixElements = this->ProcessTwoBodyMatrixElements(TmpSigma1, TmpSigma2, TmpSigma3, TmpSigma4, TmpKx1, TmpKy1, TmpKx2, TmpKy2, TmpKx3, TmpKy3, TmpKx4, TmpKy4, TmpMatrixElements);
  int* TmpLinearizedSumK = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK1 = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK2 = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK3 = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK4 = new int[TmpNbrTwoBodyMatrixElements];
  for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
    {
      TmpLinearizedSumK[i] = this->TightBindingModel->GetLinearizedMomentumIndex((TmpKx1[i] + TmpKx2[i]) % this->NbrSiteX,
										 (TmpKy1[i] + TmpKy2[i]) % this->NbrSiteY);
      TmpLinearizedK1[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx1[i], TmpKy1[i]);
      TmpLinearizedK2[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx2[i], TmpKy2[i]);
      TmpLinearizedK3[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx3[i], TmpKy3[i]);
      TmpLinearizedK4[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx4[i], TmpKy4[i]);
      TmpMatrixElements[i] *= this->TwoBodyInteractionRescalingFactor;
    }
  bool**** InternalIndicesFlags = this->TestMatrixElementsConservedDegreesOfFreedom(TmpNbrTwoBodyMatrixElements, TmpSigma1, TmpSigma2, TmpSigma3, TmpSigma4);
  
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
	  {
	    ++this->NbrInterSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY)];
	  }
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }

  int*** TmpLinearizedKInterIndices = new int** [this->NbrInterSectorSums];
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      TmpLinearizedKInterIndices[j] = new int* [this->NbrSiteX * this->NbrSiteY];
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  TmpLinearizedKInterIndices[j][i] = new int [this->NbrSiteX * this->NbrSiteY];
	}
    }
  int*** TmpLinearizedKIntraIndices = new int** [this->NbrIntraSectorSums];
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      TmpLinearizedKIntraIndices[j] = new int* [this->NbrSiteX * this->NbrSiteY];
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  TmpLinearizedKIntraIndices[j][i] = new int [this->NbrSiteX * this->NbrSiteY];
	}
    }
  
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)    
	  {
	    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
	    TmpLinearizedKInterIndices[TmpSum][this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1)][this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2)] = this->NbrInterSectorIndicesPerSum[TmpSum];
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
										     (ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    TmpLinearizedKIntraIndices[TmpSum][this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1)][this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2)] = this->NbrIntraSectorIndicesPerSum[TmpSum];
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      

      double TmpKx1 = 0.0;
      double TmpKx2 = 0.0;
      double TmpKx3 = 0.0;
      double TmpKx4 = 0.0;
      double TmpKy1 = 0.0;
      double TmpKy2 = 0.0;
      double TmpKy3 = 0.0;
      double TmpKy4 = 0.0;
      
      this->InteractionFactorsSigma = new Complex***** [this->NbrInternalIndices];
      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
	{
	  this->InteractionFactorsSigma[sigma3] = new Complex****  [this->NbrInternalIndices];
	  for (int sigma4 = sigma3; sigma4 < this->NbrInternalIndices; ++sigma4)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma4] = new Complex***[this->NbrInternalIndices];
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  this->InteractionFactorsSigma[sigma3][sigma4][sigma1] = new Complex**[this->NbrInternalIndices];
		  for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)		
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = 0;
		    }
		}
	    }
	}
      
      if (this->AdditionalSpinFlag == false)
	{
	  // spinless case
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		{
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      int Tmp = this->NbrIntraSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j];
		      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j] = new Complex [Tmp];
		    }
		}
	    }
	  
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
		      for (int j = 0; j < this->NbrInterSectorSums; ++j)
			{
			  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
			}
		    }
		}
	    }

	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		      for (int j = 0; j < this->NbrInterSectorSums; ++j)
			{
			  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
			  
			}
		    }
		}
	    }
	  
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			{
			  if (InternalIndicesFlags[sigma3][sigma4][sigma1][sigma2] == true)
			    {
			      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
			      for (int j = 0; j < this->NbrInterSectorSums; ++j)
				{
				  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
				}
			    }
			  else
			    {
			      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = 0;
			    }
			}
		    }
		}
	    }
	  
	  for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
	    {
	      int TmpIndex = 0;
	      int TmpMaxIndexFactor = 0;
	      int TmpSumK = TmpLinearizedSumK[i];
	      int Sigma1 = TmpSigma1[i];
	      int Sigma2 = TmpSigma2[i];
	      int Sigma3 = TmpSigma3[i];
	      int Sigma4 = TmpSigma4[i];
	      int K1 = TmpLinearizedK1[i];
	      int K2 = TmpLinearizedK2[i];
	      int K3 = TmpLinearizedK3[i];
	      int K4 = TmpLinearizedK4[i];
	      double TmpSign = 1.0;
	      if (Sigma2 < Sigma1)
		{
		  int Tmp = Sigma2;
		  Sigma2 = Sigma1;
		  Sigma1 = Tmp;
		  Tmp = K2;
		  K2 = K1;
		  K1 = Tmp;
		  TmpSign *= -1.0;
		}
	      if (Sigma4 < Sigma3)
		{
		  int Tmp = Sigma4;
		  Sigma4 = Sigma3;
		  Sigma3 = Tmp;
		  Tmp = K4;
		  K4 = K3;
		  K3 = Tmp;
		  TmpSign *= -1.0;
		}
	      if (Sigma1 == Sigma2)
		{
		  if (K1 > K2)
		    {
		      int Tmp = K2;
		      K2 = K1;
		      K1 = Tmp;
		      TmpSign *= -1.0;
		    }
		  if (K1 == K2)
		    {
		      TmpSign = 0.0;
		    }
		  else
		    {
		      if (Sigma3 == Sigma4)
			{
			  if (K3 > K4)
			    {
			      int Tmp = K4;
			      K4 = K3;
			      K3 = Tmp;
			      TmpSign *= -1.0;
			    }
			  if (K3 == K4)
			    {
			      TmpSign = 0.0;
			    }
			  else
			    {
			      TmpIndex = ((this->NbrIntraSectorIndicesPerSum[TmpSumK] * TmpLinearizedKIntraIndices[TmpSumK][K3][K4])
					  + TmpLinearizedKIntraIndices[TmpSumK][K1][K2]);
			    }
			}
		      else
			{
			  TmpIndex = ((this->NbrIntraSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
				      + TmpLinearizedKIntraIndices[TmpSumK][K1][K2]);
			}
		    }
		}
	      else
		{
		  if (Sigma3 == Sigma4)
		    {
		      if (K3 > K4)
			{
			  int Tmp = K4;
			  K4 = K3;
			  K3 = Tmp;
			  TmpSign *= -1.0;
			}
		      if (K3 == K4)
			{
			  TmpSign = 0.0;
			}
		      else
			{
			  TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKIntraIndices[TmpSumK][K3][K4])
				      + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
			}
		    }
		  else
		    {
		      TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
				  + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
		    }
		}
	      if (TmpSign != 0.0)
		{
		  this->InteractionFactorsSigma[Sigma1][Sigma2][Sigma3][Sigma4][TmpSumK][TmpIndex] += TmpSign * TmpMatrixElements[i];
		  ++TotalNbrInteractionFactors;
		}
	    }
	}
      else
	{
	  // spinful case
	  int ReducedNbrInternalIndices = this->NbrInternalIndices / 2;
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      for (int sigma4 = sigma3; sigma4 < this->NbrInternalIndices; ++sigma4)
			{
			  if ((InternalIndicesFlags[sigma3 % ReducedNbrInternalIndices][sigma4 % ReducedNbrInternalIndices][sigma1 % ReducedNbrInternalIndices][sigma2 % ReducedNbrInternalIndices] == true) &&
			      (this->TestSpinfulValleyConservation(sigma1, sigma2, sigma3, sigma4) == true))
			    //			      ((((sigma1 & 2) == (sigma3 & 2)) && ((sigma2 & 2) == (sigma4 & 2)))))
			    {
			      if (sigma3 == sigma4)
				{
				  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
				}
			      else
				{
				  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
				}
			      for (int j = 0; j < this->NbrInterSectorSums; ++j)
				{
				  int Tmp;
				  if (sigma3 == sigma4)
				    {
				      Tmp = this->NbrIntraSectorIndicesPerSum[j];
				    }
				  else
				    {
				      Tmp = this->NbrInterSectorIndicesPerSum[j];
				    }
				  if (sigma1 == sigma2)
				    {
				      Tmp *= this->NbrIntraSectorIndicesPerSum[j];
				    }
				  else
				    {
				      Tmp *= this->NbrInterSectorIndicesPerSum[j];
				    }
				  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j] = new Complex [Tmp];
				  Complex* TmpInteractionArray = this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j];
				}
			    }
			}
		    }
		}
	    }

	  for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
	    {
	      int TmpIndex = 0;
	      int TmpMaxIndexFactor = 0;
	      int TmpSumK = TmpLinearizedSumK[i];
	      int Sigma1 = TmpSigma1[i];
	      int Sigma2 = TmpSigma2[i];
	      int Sigma3 = TmpSigma3[i];
	      int Sigma4 = TmpSigma4[i];
	      int K1 = TmpLinearizedK1[i];
	      int K2 = TmpLinearizedK2[i];
	      int K3 = TmpLinearizedK3[i];
	      int K4 = TmpLinearizedK4[i];
	      double TmpSign = 1.0;
	      if (Sigma2 < Sigma1)
		{
		  int Tmp = Sigma2;
		  Sigma2 = Sigma1;
		  Sigma1 = Tmp;
		  Tmp = K2;
		  K2 = K1;
		  K1 = Tmp;
		  TmpSign *= -1.0;
		}
	      if (Sigma4 < Sigma3)
		{
		  int Tmp = Sigma4;
		  Sigma4 = Sigma3;
		  Sigma3 = Tmp;
		  Tmp = K4;
		  K4 = K3;
		  K3 = Tmp;
		  TmpSign *= -1.0;
		}
	      if (Sigma1 == Sigma2)
		{
		  if (K1 > K2)
		    {
		      int Tmp = K2;
		      K2 = K1;
		      K1 = Tmp;
		      TmpSign *= -1.0;
		    }
		  if (K1 == K2)
		    {
		      TmpSign = 0.0;
		    }
		  else
		    {
		      if (Sigma3 == Sigma4)
			{
			  if (K3 > K4)
			    {
			      int Tmp = K4;
			      K4 = K3;
			      K3 = Tmp;
			      TmpSign *= -1.0;
			    }
			  if (K3 == K4)
			    {
			      TmpSign = 0.0;
			    }
			  else
			    {
			      TmpIndex = ((this->NbrIntraSectorIndicesPerSum[TmpSumK] * TmpLinearizedKIntraIndices[TmpSumK][K3][K4])
					  + TmpLinearizedKIntraIndices[TmpSumK][K1][K2]);
			    }
			}
		      else
			{
			  TmpIndex = ((this->NbrIntraSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
				      + TmpLinearizedKIntraIndices[TmpSumK][K1][K2]);
			}
		    }
		}
	      else
		{
		  if (Sigma3 == Sigma4)
		    {
		      if (K3 > K4)
			{
			  int Tmp = K4;
			  K4 = K3;
			  K3 = Tmp;
			  TmpSign *= -1.0;
			}
		      if (K3 == K4)
			{
			  TmpSign = 0.0;
			}
		      else
			{
			  TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKIntraIndices[TmpSumK][K3][K4])
				      + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
			}
		    }
		  else
		    {
		      TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
				  + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
		    }
		}
	      if (TmpSign != 0.0)
		{
		  this->InteractionFactorsSigma[Sigma1][Sigma2][Sigma3][Sigma4][TmpSumK][TmpIndex] += TmpSign * TmpMatrixElements[i];
		  this->InteractionFactorsSigma[Sigma1 + ReducedNbrInternalIndices][Sigma2 + ReducedNbrInternalIndices][Sigma3 + ReducedNbrInternalIndices][Sigma4 + ReducedNbrInternalIndices][TmpSumK][TmpIndex] += TmpSign * TmpMatrixElements[i];
		  TotalNbrInteractionFactors += 2l;
		}
	    }	  	
	  for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
	    {
	      int TmpIndex = 0;
	      int TmpMaxIndexFactor = 0;
	      int TmpSumK = TmpLinearizedSumK[i];
	      int Sigma1 = TmpSigma1[i];
	      int Sigma2 = TmpSigma2[i];
	      int Sigma3 = TmpSigma3[i];
	      int Sigma4 = TmpSigma4[i];
	      int K1 = TmpLinearizedK1[i];
	      int K2 = TmpLinearizedK2[i];
	      int K3 = TmpLinearizedK3[i];
	      int K4 = TmpLinearizedK4[i];	      
	      double TmpSign = 1.0;
	      TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K4][K3])
	       		  + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
	      TmpSign = -1.0;
	      this->InteractionFactorsSigma[Sigma1][Sigma2 + ReducedNbrInternalIndices][Sigma4][Sigma3 + ReducedNbrInternalIndices][TmpSumK][TmpIndex] += TmpSign * TmpMatrixElements[i];
	      TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
	       		  + TmpLinearizedKInterIndices[TmpSumK][K2][K1]);
	      TmpSign = -1.0;
	      this->InteractionFactorsSigma[Sigma2][Sigma1 + ReducedNbrInternalIndices][Sigma3][Sigma4 + ReducedNbrInternalIndices][TmpSumK][TmpIndex] += TmpSign * TmpMatrixElements[i];
	      TotalNbrInteractionFactors += 2l;
	    }
	}
    }
  else
    {      
      // bosonic interaction
    }

  this->FreeMatrixElementsConservedDegreesOfFreedom(InternalIndicesFlags);
  delete[] TmpSigma1;
  delete[] TmpSigma2;
  delete[] TmpSigma3;
  delete[] TmpSigma4;
  delete[] TmpKx1;
  delete[] TmpKx2;
  delete[] TmpKx3;
  delete[] TmpKx4;
  delete[] TmpKy1;
  delete[] TmpKy2;
  delete[] TmpKy3;
  delete[] TmpKy4;
  delete[] TmpMatrixElements;
  delete[] TmpLinearizedSumK;
  delete[] TmpLinearizedK1;
  delete[] TmpLinearizedK2;
  delete[] TmpLinearizedK3;
  delete[] TmpLinearizedK4;
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  delete[] TmpLinearizedKInterIndices[j][i];
	}
      delete[] TmpLinearizedKInterIndices[j];
    }
  delete[] TmpLinearizedKInterIndices;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  delete[] TmpLinearizedKIntraIndices[j][i];
	}
      delete[] TmpLinearizedKIntraIndices[j];
    }
  delete[] TmpLinearizedKIntraIndices;
   cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// evaluate all one-body factors
//   

void ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::EvaluateOneBodyInteractionFactors()
{
  this->OneBodyInteractionFactorsSigma = new Complex**[this->NbrInternalIndices];
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      this->OneBodyInteractionFactorsSigma[sigma1] = new Complex*[this->NbrInternalIndices];
      for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
	{
	  this->OneBodyInteractionFactorsSigma[sigma1][sigma2] = 0;
	}
    }
      
  if (this->FlatBand == false)
    { 
      bool** TmpOneBodyFlags = new bool*[this->NbrInternalIndices];
      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  TmpOneBodyFlags[sigma1] = new bool[this->NbrInternalIndices];
	  for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      TmpOneBodyFlags[sigma1][sigma2] = false;
	    }
	}
      
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      HermitianMatrix TmpBlochHamiltonian(this->TightBindingModel->ComputeBlochHamiltonian(kx, ky));
	      Complex Tmp;
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
		    {
		      TmpBlochHamiltonian.GetMatrixElement(sigma1, sigma2, Tmp);
		      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
			{
			  if (TmpOneBodyFlags[sigma1][sigma2] == false)
			    {
			      this->OneBodyInteractionFactorsSigma[sigma1][sigma2] = new Complex [this->TightBindingModel->GetNbrStatePerBand()];			      
			      TmpOneBodyFlags[sigma1][sigma2] = true;
			    }
			}
		    }
		}
	    }
	}
      
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      HermitianMatrix TmpBlochHamiltonian(this->TightBindingModel->ComputeBlochHamiltonian(kx, ky));
	      Complex Tmp;
	      int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
		    {
		      if (TmpOneBodyFlags[sigma1][sigma2] == true)
			{
			  TmpBlochHamiltonian.GetMatrixElement(sigma1, sigma2, Tmp);
			  this->OneBodyInteractionFactorsSigma[sigma1][sigma2][Index] = Tmp;
			}
		    }
		}
	    }
	}
    }
}

// evaluate all one-body factors
// process the three-body matrix elements from the ascii file
//
// arraySigma1 = reference on the array containing the indices of the internal degree freedom for the operator 1
// arraySigma2 = reference on the array containing the indices of the internal degree freedom for the operator 2
// arraySigma3 = reference on the array containing the indices of the internal degree freedom for the operator 3
// arraySigma4 = reference on the array containing the indices of the internal degree freedom for the operator 4
// arraySigma5 = reference on the array containing the indices of the internal degree freedom for the operator 5
// arraySigma6 = reference on the array containing the indices of the internal degree freedom for the operator 6
// arrayKx1 = reference on the array containing the momentum along x for the operator 1
// arrayKy1 = reference on the array containing the momentum along y for the operator 1
// arrayKx2 = reference on the array containing the momentum along x for the operator 2
// arrayKy2 = reference on the array containing the momentum along y for the operator 2
// arrayKx3 = reference on the array containing the momentum along x for the operator 3
// arrayKy3 = reference on the array containing the momentum along y for the operator 3
// arrayKx4 = reference on the array containing the momentum along x for the operator 4
// arrayKy4 = reference on the array containing the momentum along y for the operator 4
// arrayKx5 = reference on the array containing the momentum along x for the operator 5
// arrayKy5 = reference on the array containing the momentum along y for the operator 5
// arrayKx6 = reference on the array containing the momentum along x for the operator 6
// arrayKy6 = reference on the array containing the momentum along y for the operator 6
// arrayMatrixElements = reference on the array containing the matrix elements
// return value = number of entries

int ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::ProcessThreeBodyMatrixElements(int*& arraySigma1, int*& arraySigma2, int*& arraySigma3, int*& arraySigma4, int*& arraySigma5, int*& arraySigma6,
												    int*& arrayKx1, int*& arrayKy1, int*& arrayKx2, int*& arrayKy2, int*& arrayKx3, int*& arrayKy3,
												    int*& arrayKx4, int*& arrayKy4, int*& arrayKx5, int*& arrayKy5, int*& arrayKx6, int*& arrayKy6, Complex*& arrayMatrixElements)
{
  MultiColumnASCIIFile TmpInteractionFile;
  if (TmpInteractionFile.Parse(this->MatrixElementsThreeBodyInteractionFile) == false)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrLines() == 0)
    {
      cout << this->MatrixElementsThreeBodyInteractionFile << " is an empty file" << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrColumns() < 13)
    {
      cout << this->MatrixElementsThreeBodyInteractionFile << " has a wrong number of column (has "
	   << TmpInteractionFile.GetNbrColumns() << ", should be at least 13)" << endl;
      exit(0);
    }
  int TmpNbrThreeBodyMatrixElements = TmpInteractionFile.GetNbrLines();
  cout << "nbr of two body matrix elements in " << this->MatrixElementsThreeBodyInteractionFile << " = " << TmpNbrThreeBodyMatrixElements << endl;

  arraySigma1 = new int [TmpNbrThreeBodyMatrixElements];
  arraySigma2 = new int [TmpNbrThreeBodyMatrixElements];
  arraySigma3 = new int [TmpNbrThreeBodyMatrixElements];
  arraySigma4 = new int [TmpNbrThreeBodyMatrixElements];
  arraySigma5 = new int [TmpNbrThreeBodyMatrixElements];
  arraySigma6 = new int [TmpNbrThreeBodyMatrixElements];
  arrayKx1 = TmpInteractionFile.GetAsIntegerArray(0);
  arrayKx2 = TmpInteractionFile.GetAsIntegerArray(2);
  arrayKx3 = TmpInteractionFile.GetAsIntegerArray(4);
  arrayKx4 = TmpInteractionFile.GetAsIntegerArray(6);
  arrayKx5 = TmpInteractionFile.GetAsIntegerArray(8);
  arrayKx6 = TmpInteractionFile.GetAsIntegerArray(10);
  arrayKy1 = TmpInteractionFile.GetAsIntegerArray(1);
  arrayKy2 = TmpInteractionFile.GetAsIntegerArray(3);
  arrayKy3 = TmpInteractionFile.GetAsIntegerArray(5);
  arrayKy4 = TmpInteractionFile.GetAsIntegerArray(7);
  arrayKy5 = TmpInteractionFile.GetAsIntegerArray(9);
  arrayKy6 = TmpInteractionFile.GetAsIntegerArray(11);
  arrayMatrixElements = TmpInteractionFile.GetAsComplexArray(12);
  for (int i = 0; i < TmpNbrThreeBodyMatrixElements; ++i)
    {
      arraySigma1[i] = 0;
      arraySigma2[i] = 0;
      arraySigma3[i] = 0;
      arraySigma4[i] = 0;
      arraySigma5[i] = 0;
      arraySigma6[i] = 0;
    }
  if (arrayMatrixElements == 0)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  return TmpNbrThreeBodyMatrixElements;
}

// process the matrix elements from the ascii file
//
// arraySigma1 = reference on the array containing the indices of the internal degree freedom for the operator 1
// arraySigma2 = reference on the array containing the indices of the internal degree freedom for the operator 2
// arraySigma3 = reference on the array containing the indices of the internal degree freedom for the operator 3
// arraySigma4 = reference on the array containing the indices of the internal degree freedom for the operator 4
// arrayKx1 = reference on the array containing the momentum along x for the operator 1
// arrayKy1 = reference on the array containing the momentum along y for the operator 1
// arrayKx2 = reference on the array containing the momentum along x for the operator 2
// arrayKy2 = reference on the array containing the momentum along y for the operator 2
// arrayKx3 = reference on the array containing the momentum along x for the operator 3
// arrayKy3 = reference on the array containing the momentum along y for the operator 3
// arrayKx4 = reference on the array containing the momentum along x for the operator 4
// arrayKy4 = reference on the array containing the momentum along y for the operator 4
// arrayMatrixElements = reference on the array containing the matrix elements
// return value = number of entries

int ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian::ProcessTwoBodyMatrixElements(int*& arraySigma1, int*& arraySigma2, int*& arraySigma3, int*& arraySigma4, int*& arrayKx1, int*& arrayKy1, int*& arrayKx2, int*& arrayKy2, int*& arrayKx3, int*& arrayKy3, int*& arrayKx4, int*& arrayKy4, Complex*& arrayMatrixElements)
{
  MultiColumnASCIIFile TmpInteractionFile;
  if (TmpInteractionFile.Parse(this->MatrixElementsTwoBodyInteractionFile) == false)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrLines() == 0)
    {
      cout << this->MatrixElementsTwoBodyInteractionFile << " is an empty file" << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrColumns() < 9)
    {
      cout << this->MatrixElementsTwoBodyInteractionFile << " has a wrong number of column (has "
	   << TmpInteractionFile.GetNbrColumns() << ", should be at least 9)" << endl;
      exit(0);
    }
  int TmpNbrTwoBodyMatrixElements = TmpInteractionFile.GetNbrLines();
  cout << "nbr of two body matrix elements in " << this->MatrixElementsTwoBodyInteractionFile << " = " << TmpNbrTwoBodyMatrixElements << endl;

  arraySigma1 = new int [TmpNbrTwoBodyMatrixElements];
  arraySigma2 = new int [TmpNbrTwoBodyMatrixElements];
  arraySigma3 = new int [TmpNbrTwoBodyMatrixElements];
  arraySigma4 = new int [TmpNbrTwoBodyMatrixElements];
  arrayKx1 = TmpInteractionFile.GetAsIntegerArray(0);
  arrayKx2 = TmpInteractionFile.GetAsIntegerArray(2);
  arrayKx3 = TmpInteractionFile.GetAsIntegerArray(4);
  arrayKx4 = TmpInteractionFile.GetAsIntegerArray(6);
  arrayKy1 = TmpInteractionFile.GetAsIntegerArray(1);
  arrayKy2 = TmpInteractionFile.GetAsIntegerArray(3);
  arrayKy3 = TmpInteractionFile.GetAsIntegerArray(5);
  arrayKy4 = TmpInteractionFile.GetAsIntegerArray(7);
  arrayMatrixElements = TmpInteractionFile.GetAsComplexArray(8);
  for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
    {
      arraySigma1[i] = 0;
      arraySigma2[i] = 0;
      arraySigma3[i] = 0;
      arraySigma4[i] = 0;
    }
  if (arrayMatrixElements == 0)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  return TmpNbrTwoBodyMatrixElements;
}

