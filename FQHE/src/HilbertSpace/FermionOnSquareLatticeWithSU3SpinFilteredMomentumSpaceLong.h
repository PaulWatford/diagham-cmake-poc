////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(3) spin         //
//                  for more than 21 orbitals in momentum space               //
//                          in a filtered momentum space                      //
//                                                                            //
//                        last modification : 14/12/2023                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU3SPINFILTEREDMOMENTUMSPACELONG_H
#define FERMIONONSQUARELATTICEWITHSU3SPINFILTEREDMOMENTUMSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong.h"

#include <iostream>



class FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong : public FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong
{

 protected:

  // orbital filering mask (i.e. state should be 0 once the mask is applied)
  ULONGLONG OrbitalFilteringMask;
  
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
  FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, char* allowedOrbitalsFileName, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong(const FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong& operator = (const FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // filter Hilbert to remove forbidden orbitals
  //
  // allowedOrbitalsFileName = ascii file providing the orbitals that are allowed
  virtual void FilterHilbertSpace(char* allowedOrbitalsFileName);

};


#endif


