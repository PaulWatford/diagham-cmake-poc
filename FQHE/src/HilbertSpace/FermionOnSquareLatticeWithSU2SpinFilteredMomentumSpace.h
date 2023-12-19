////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on square lattice with spin            //
//                        in filtered momentum space                          //
// using a simplified version more in line with SU(3), SU(6) or SU(12) codes  //
//                                                                            //
//                        last modification : 17/12/2023                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU2SPINFILTEREDMOMENTUMSPACE_H
#define FERMIONONSQUARELATTICEWITHSU2SPINFILTEREDMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include <iostream>



class FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace : public FermionOnSquareLatticeWithSU2SpinMomentumSpace
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:

  // orbital filering mask (i.e. state should be 0 once the mask is applied)
  unsigned long OrbitalFilteringMask;

public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, char* allowedOrbitalsFileName, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // nbrSpinUp = number of particles with spin up
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, char* allowedOrbitalsFileName, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace(const FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace& operator = (const FermionOnSquareLatticeWithSU2SpinFilteredMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

protected:

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateFilteredStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateFilteredHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);
  
  // parse the ascii file providing the orbitals that are allowed
  //
  // allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
  virtual void ParseOrbitalFile(char* allowedOrbitalsFileName);

  // filter Hilbert to remove forbidden orbitals
  //
  // allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
  virtual void FilterHilbertSpace(char* allowedOrbitalsFileName);

};


#endif


