////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                 class of fermions on sphere with N flavors and             //
//                        Lz<->-Lz symmetry (aka inversion)                   //
//                                                                            //
//                           class author: Sahana Das                         //
//                                                                            //
//                        last modification : 27/05/2026                      //
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



#ifndef FERMIONONSPHEREWITHNFLAVORLZSYMMETRY_H
#define FERMIONONSPHEREWITHNFLAVORLZSYMMETRY_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithNFlavor.h"

#include <iostream>
#include <cstring>
#include <cmath>

using std::cout;
using std::endl;


class FermionOnSphereWithNFlavorLzSymmetry : public FermionOnSphereWithNFlavor
{
  
public:

  // Keep a 6-arg interface (legacy style) so your program/manager can stay unchanged.
  // nbrFermions and memory are currently unused here because the base class
  // computes N from nbrParticlesPerFlavor and ParticleOnSphere has no such ctor.
  FermionOnSphereWithNFlavorLzSymmetry(int nbrFermions,
                                       int totalLz,
                                       int lzMax,
                                       int nbrFlavors,
                                       int* nbrParticlesPerFlavor,
                                       unsigned long memory);

  virtual ~FermionOnSphereWithNFlavorLzSymmetry();

  virtual AbstractHilbertSpace* Clone();
  virtual std::ostream& PrintState(std::ostream& os, int state);
  virtual int GetParticleStatistic();

  // This matches the style used across DiagHam (optional, but harmless).
  virtual int GetHilbertSpaceAdditionalSymmetry();

protected:
  unsigned long Memory; // stored for compatibility; not required by current implementation

private:
  FermionOnSphereWithNFlavorLzSymmetry();
};

#endif
