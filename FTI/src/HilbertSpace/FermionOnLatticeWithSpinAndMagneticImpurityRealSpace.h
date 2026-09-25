////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//                  in real space with one magnetic impurity                  //
//                                                                            //
//                        last modification : 19/08/2022                      //
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


#ifndef FERMIONONLATTICEWITHSPINANDMAGNETICIMPURITYREALSPACE_H
#define FERMIONONLATTICEWITHSPINANDMAGNETICIMPURITYREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"

#include <iostream>



class FermionOnLatticeWithSpinAndMagneticImpurityRealSpace : public FermionOnLatticeWithSpinRealSpace
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;
  friend class FermionOnLatticeWithSpinRealSpaceAnd2DTranslation;
  friend class FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation;

 protected:


 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinAndMagneticImpurityRealSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinAndMagneticImpurityRealSpace (int nbrFermions, int nbrSite, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the total spin value
  // nbrSite = number of sites in the x direction
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinAndMagneticImpurityRealSpace (int nbrFermions, int totalSpin, int nbrSite, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinAndMagneticImpurityRealSpace(const FermionOnLatticeWithSpinAndMagneticImpurityRealSpace& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinAndMagneticImpurityRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinAndMagneticImpurityRealSpace& operator = (const FermionOnLatticeWithSpinAndMagneticImpurityRealSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real space
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);

};


#endif


