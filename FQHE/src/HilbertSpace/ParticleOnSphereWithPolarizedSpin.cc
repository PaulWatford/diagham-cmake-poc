////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of particle on sphere with spin                  //
//                 using internally a fully polarized Hilbert space           //
//                                                                            //
//                        last modification : 03/01/2024                      //
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
#include "HilbertSpace/ParticleOnSphereWithPolarizedSpin.h"


// default constructor
//

ParticleOnSphereWithPolarizedSpin::ParticleOnSphereWithPolarizedSpin()
{
  this->PolarizedHilbertSpace = 0;
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0l;
}

// constructor
//
// polarizedHilbertSpace = pointer to the polarized Hilbert space

ParticleOnSphereWithPolarizedSpin::ParticleOnSphereWithPolarizedSpin(ParticleOnSphere* polarizedHilbertSpace)
{
  this->PolarizedHilbertSpace = (ParticleOnSphere*) polarizedHilbertSpace->Clone();
  this->HilbertSpaceDimension = this->PolarizedHilbertSpace->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->PolarizedHilbertSpace->GetLargeHilbertSpaceDimension();
}
  
// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

ParticleOnSphereWithPolarizedSpin::ParticleOnSphereWithPolarizedSpin(const ParticleOnSphereWithPolarizedSpin& fermions)
{
  this->PolarizedHilbertSpace = (ParticleOnSphere*) (fermions.PolarizedHilbertSpace->Clone());
  this->HilbertSpaceDimension = this->PolarizedHilbertSpace->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->PolarizedHilbertSpace->GetLargeHilbertSpaceDimension();
}

// destructor
//

ParticleOnSphereWithPolarizedSpin::~ParticleOnSphereWithPolarizedSpin ()
{
  if (this->PolarizedHilbertSpace != 0)
    {
      delete this->PolarizedHilbertSpace;
    }
}

// assignement (without duplicating data)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

ParticleOnSphereWithPolarizedSpin& ParticleOnSphereWithPolarizedSpin::operator = (const ParticleOnSphereWithPolarizedSpin& fermions)
{
  this->PolarizedHilbertSpace = (ParticleOnSphere*) (fermions.PolarizedHilbertSpace->Clone());
  this->HilbertSpaceDimension = this->PolarizedHilbertSpace->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->PolarizedHilbertSpace->GetLargeHilbertSpaceDimension();
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithPolarizedSpin::Clone()
{
  return new ParticleOnSphereWithPolarizedSpin(*this);
}
