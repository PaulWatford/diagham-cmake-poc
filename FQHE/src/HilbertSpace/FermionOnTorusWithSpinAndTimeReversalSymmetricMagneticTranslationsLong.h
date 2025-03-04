////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of fermion with spin on a torus with time            //
//         reversal symmetry, taking into account magnetic translations       //
//                                                                            //
//                        last modification : 04/03/2025                      //
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


#ifndef FERMIONONTORUSWITHSPINANDTIMEREVERSALSYMMETRICMAGNETICTRANSLATIONSLONG_H
#define FERMIONONTORUSWITHSPINANDTIMEREVERSALSYMMETRICMAGNETICTRANSLATIONSLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslationsLong.h"

using std::cout;
using std::endl;
using std::dec;
using std::hex;




class FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong :  public FermionOnTorusWithSpinAndMagneticTranslationsLong
{  
   public:
     
  // basic constructor
  // 
  // nbrFermions= number of fermions
  // totalSpin = twice the total spin value
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)  
  FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong (int nbrFermions, int totalSpin, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong(const FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong& fermions);

  // destructor
  //
  ~FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong& operator = (const FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();
  
  protected:
    
    
  // evaluate Hilbert space dimension using recursive algorithm
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalMomentum = momentum total value
  // totalSpin = number of particles with spin up
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalMomentum, int totalSpin);
      
  // evaluate Hilbert space dimension without the translation symmmetry along x
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalMomentum = momentum total value
  // totalSpinUp = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = Hilbert space dimension
  virtual long RawGenerateStates(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp, long pos);
  
 
  // get the total spin
  //
  //return value: total spin of the Hilbert space
  virtual int GetTotalSpin();

};

// get the total spin
//
//return value: total spin of the Hilbert space
inline int FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslationsLong::GetTotalSpin()
{
 return this->TotalSpin; 
}
#endif
