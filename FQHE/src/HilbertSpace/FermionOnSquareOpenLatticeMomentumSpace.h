////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on square lattice                 //
//            with open boundary conditions along the y direction             //
//                   and in momentum space for the x direction                //
//                                                                            //
//                        last modification : 07/08/2026                      //
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


#ifndef FERMIONONOPENSQUARELATTICEMOMENTUMSPACE_H
#define FERMIONONOPENSQUARELATTICEMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include <iostream>



class FermionOnSquareOpenLatticeMomentumSpace : public FermionOnSquareLatticeMomentumSpace
{

  friend class FermionOnSquareLatticeMomentumSpace;
  friend class FermionOnSquareLatticeWithSU2SpinMomentumSpace;
  friend class FermionOnSquareLatticeWithSU3SpinMomentumSpace;
  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;
  friend class FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong;
  friend class FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong;
  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong;
  
 protected:


 public:

  // default constructor
  // 
  FermionOnSquareOpenLatticeMomentumSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // outputDirectory = if non-zero, the constructor looks for a previously saved Hilbert space and if not avaliable, will save the current one after generation
  // memory = amount of memory granted for precalculations
  FermionOnSquareOpenLatticeMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, char* outputDirectory = 0, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareOpenLatticeMomentumSpace(const FermionOnSquareOpenLatticeMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareOpenLatticeMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareOpenLatticeMomentumSpace& operator = (const FermionOnSquareOpenLatticeMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // provide the default name of the file for Hilbert space storage 
  //
  // return value = pointer to file name (0 if no default file name exists) 
  virtual char* GetDefaultHilbertSpaceFileName();

protected:

  // core part of the Hilbert space generation
  //
  void GenerateCoreHilbertSpace();
    
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current position along y for a single particle
  // currentTotalKx = current total momentum along x
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current position along y for a single particle
  // currentTotalKx = current total momentum along x
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, long pos);

};


#endif


