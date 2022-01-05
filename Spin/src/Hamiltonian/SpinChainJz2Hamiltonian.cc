////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of spin chain hamiltonian                      //
//                    with an additional on-site Jz^2 term                    //
//                                                                            //
//                        last modification : 07/12/2021                      //
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


#include "Hamiltonian/SpinChainJz2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = array containing coupling constants between spins
// jz2 = array containing the amplitude of the on-site Jz^2 term
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainJz2Hamiltonian::SpinChainJz2Hamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz2, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->Hz = 0;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->J = new double [this->NbrSpin - 1];
      this->Jz = new double [this->NbrSpin - 1];
      this->HalfJ = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = j[i];
	  this->HalfJ[i] = j[i] * 0.5;
	}
    }
  else
    {
      this->J = new double [this->NbrSpin];
      this->Jz = new double [this->NbrSpin];
      this->HalfJ = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = j[i];
	  this->HalfJ[i] = j[i] * 0.5;
	}
    }
  this->Jz2 = new double [this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; i++)
    {
      this->Jz2[i] = jz2[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = array containing coupling constants between spins along x and z
// jz = array containing coupling constants between spins along z
// jz2 = array containing the amplitude of the on-site Jz^2 term
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainJz2Hamiltonian::SpinChainJz2Hamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* jz2, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->Hz = 0;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->J = new double [this->NbrSpin - 1];
      this->Jz = new double [this->NbrSpin - 1];
      this->HalfJ = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	}
    }
  else
    {
      cout << "check " << endl;
      this->J = new double [this->NbrSpin];
      this->Jz = new double [this->NbrSpin];
      this->HalfJ = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	}
    }
  this->Jz2 = new double [this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; i++)
    {
      this->Jz2[i] = jz2[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// constructor with a generic magnetic field
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = array containing coupling constants between spins along x and z
// jz = array containing coupling constants between spins along z
// hz = array containing the amplitude of the Zeeman term along z
// jz2 = array containing the amplitude of the on-site Jz^2 term
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainJz2Hamiltonian::SpinChainJz2Hamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* hz, double* jz2, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->J = new double [this->NbrSpin - 1];
      this->Jz = new double [this->NbrSpin - 1];
      this->HalfJ = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	}
    }
  else
    {
      this->J = new double [this->NbrSpin];
      this->Jz = new double [this->NbrSpin];
      this->HalfJ = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	}
   }
  this->Hz = new double [this->NbrSpin];
  this->Jz2 = new double [this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; i++)
    {
      this->Hz[i] = hz[i];
      this->Jz2[i] = jz2[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SpinChainJz2Hamiltonian::~SpinChainJz2Hamiltonian() 
{
  delete[] this->Jz2;
}

// evaluate diagonal matrix elements
// 

void SpinChainJz2Hamiltonian::EvaluateDiagonalMatrixElements()
{
  this->SpinChainHamiltonian::EvaluateDiagonalMatrixElements();
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Coefficient;
  // Jz^2 part
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < this->NbrSpin; j++)
	{
	  this->Chain->Szi(j, i, Coefficient);
	  this->SzSzContributions[i] += this->Jz2[j] * (Coefficient * Coefficient);
	}
    }
}

