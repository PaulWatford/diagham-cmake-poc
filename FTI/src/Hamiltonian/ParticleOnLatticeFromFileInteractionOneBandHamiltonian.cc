////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//        class of a two body interaction projected onto a sinle band         //
//         from an ASCII file providing the two body matrix elements          //
//                                                                            //
//                        last modification : 27/12/2023                      //
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
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionOneBandHamiltonian.h"
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

ParticleOnLatticeFromFileInteractionOneBandHamiltonian::ParticleOnLatticeFromFileInteractionOneBandHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// matrixElementsInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// interactionRescalingFactor = global rescaling factor for the two-body interaction term
// spinFlag = include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
// onlyPreserveKx = assume that the system is periodic only along x, treating ky as a coordinate along y
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeFromFileInteractionOneBandHamiltonian::ParticleOnLatticeFromFileInteractionOneBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles,
													       int nbrSiteX, int nbrSiteY,
													       char* matrixElementsInteractionFile,
													       Abstract2DTightBindingModel* tightBindingModel, 
													       bool flatBandFlag, double interactionRescalingFactor, bool spinFlag, bool onlyPreserveKx,
													       AbstractArchitecture* architecture, long memory)
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
  this->InteractionRescalingFactor = interactionRescalingFactor;
  this->AdditionalSpinFlag = spinFlag;
  this->OnlyPreserveKx = onlyPreserveKx;
  if (this->AdditionalSpinFlag == true)
    {
      this->NbrInternalIndices = 2;
    }
  else
    {
      this->NbrInternalIndices = 1;
    }
  this->MatrixElementsInteractionFile = new char[strlen(matrixElementsInteractionFile) + 1];
  strcpy(this->MatrixElementsInteractionFile, matrixElementsInteractionFile);
  
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
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

ParticleOnLatticeFromFileInteractionOneBandHamiltonian::~ParticleOnLatticeFromFileInteractionOneBandHamiltonian()
{
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

int ParticleOnLatticeFromFileInteractionOneBandHamiltonian::ProcessTwoBodyMatrixElements(int*& arraySigma1, int*& arraySigma2, int*& arraySigma3, int*& arraySigma4, int*& arrayKx1, int*& arrayKy1, int*& arrayKx2, int*& arrayKy2, int*& arrayKx3, int*& arrayKy3, int*& arrayKx4, int*& arrayKy4, Complex*& arrayMatrixElements)
{
  MultiColumnASCIIFile TmpInteractionFile;
  if (TmpInteractionFile.Parse(this->MatrixElementsInteractionFile) == false)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrLines() == 0)
    {
      cout << this->MatrixElementsInteractionFile << " is an empty file" << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrColumns() < 9)
    {
      cout << this->MatrixElementsInteractionFile << " has a wrong number of column (has "
	   << TmpInteractionFile.GetNbrColumns() << ", should be at least 9)" << endl;
      exit(0);
    }
  int TmpNbrTwoBodyMatrixElements = TmpInteractionFile.GetNbrLines();
  cout << "nbr of two body matrix elements in " << this->MatrixElementsInteractionFile << " = " << TmpNbrTwoBodyMatrixElements << endl;

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

