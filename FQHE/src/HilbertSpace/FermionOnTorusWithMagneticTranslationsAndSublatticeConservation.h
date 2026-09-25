////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//           and conservation of the number of particles per sublattice       //
//                                                                            //
//                        last modification : 05/06/2022                      //
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


#ifndef FERMIONONTORUSWITHMAGNETICTRANSLATIONSANDSUBLATTICECONSERVATION_H
#define FERMIONONTORUSWITHMAGNETICTRANSLATIONSANDSUBLATTICECONSERVATION_H


#include "config.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnSphere.h"

using std::cout;
using std::endl;
using std::dec;
using std::hex;


class FermionOnTorusWithMagneticTranslationsAndSublatticeConservation :  public FermionOnTorusWithMagneticTranslations
{

protected:
  
  
  // number of particles occupying orbitals with an even momentum
  int EvenMomentumNbrParticles;
  
public:
  
  // default constructor
  // 
  FermionOnTorusWithMagneticTranslationsAndSublatticeConservation ();
  
  // basic constructor
  // 
  // nbrFermions = number of fermions 
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  // evenMomentumNbrParticles = number of particles occupying orbitals with an even momentum
  FermionOnTorusWithMagneticTranslationsAndSublatticeConservation (int nbrFermions, int maxMomentum, int xMomentum, int yMomentum, int evenMomentumNbrParticles);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithMagneticTranslationsAndSublatticeConservation(const FermionOnTorusWithMagneticTranslationsAndSublatticeConservation& fermions);

  // destructor
  //
  ~FermionOnTorusWithMagneticTranslationsAndSublatticeConservation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithMagneticTranslationsAndSublatticeConservation& operator = (const FermionOnTorusWithMagneticTranslationsAndSublatticeConservation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
 
 protected:

  // generate all states corresponding to the constraints (without taking into the canonical form) 
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentYMomentum = current value of the momentum in the y direction
  // currentEvenMomentumNbrParticles = current number of particles occupying orbitals with an even momentum
  // return value = position from which new states have to be stored
  long RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, long pos, int currentYMomentum, int currentEvenMomentumNbrParticles);

  };


#endif


