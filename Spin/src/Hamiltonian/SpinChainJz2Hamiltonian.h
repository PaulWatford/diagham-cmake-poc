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


#ifndef SPINCHAINJZ2HAMILTONIAN_H
#define SPINCHAINJZ2HAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/SpinChainHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class SpinChainJz2Hamiltonian : public SpinChainHamiltonian
{

 protected:
  
  // array containing the amplitude of the on-site Jz^2 term
  double* Jz2;

  public:

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = array containing coupling constants between spins
  // jz2 = array containing the amplitude of the on-site Jz^2 term
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainJz2Hamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz2, bool periodicBoundaryConditions = false);

  // constructor from default data
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = array containing coupling constants between spins along x and z
  // jz = array containing coupling constants between spins along z
  // jz2 = array containing the amplitude of the on-site Jz^2 term
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainJz2Hamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* jz2, bool periodicBoundaryConditions = false);

  // constructor with a generic magnetic field
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = array containing coupling constants between spins along x and z
  // jz = array containing coupling constants between spins along z
  // hz = array containing the amplitude of the Zeeman term along z
  // jz2 = array containing the amplitude of the on-site Jz^2 term
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainJz2Hamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* hz, double* jz2, bool periodicBoundaryConditions = false);

  // destructor
  //
  ~SpinChainJz2Hamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
