////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on sphere with N flavors                 //
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


#ifndef PARTICLEONSPHEREWITHNFLAVOR_H
#define PARTICLEONSPHEREWITHNFLAVOR_H

#include "config.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/ParticleOnSphere.h"

class ParticleOnSphereWithNFlavor : public ParticleOnSphere
{
 protected:
  int NbrFlavors;

 public:

  virtual ~ParticleOnSphereWithNFlavor();

  // particle statistic (fermion = 1)
  virtual int GetParticleStatistic() = 0;
  
  // flavor information
  inline int GetNbrFlavors() const { return NbrFlavors; }
  //virtual int GetNbrFlavors() const = 0;
  
  // apply Prod_i a^+_mi Prod_i a_ni
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply sum_sigma a^+_{m sigma} a_{m sigma}
  virtual double AdA (int index, int m);
  virtual double AdA (long index, int m);

  // apply a^+_{m sigma} a_{m sigma}
  virtual double AdsigmaAsigma (int index, int m, int sigma) = 0;

  // apply a^+_{m sigma} a_{n sigma'}
  virtual int AdsigmaAsigma (int index, int m, int n, int sigma_m, int sigma_n, double& coefficient) = 0;

  // apply a_{n1 sigma1} a_{n2 sigma2}
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2) = 0;

  // apply a^+_{m1 sigma1} a^+_{m2 sigma2}
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient) = 0;

  // Wave function evaluation (identical to SU4)
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis);
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int firstComponent, int nbrComponent);

  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int nextCoordinates);
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int nextCoordinates, int firstComponent, int nbrComponent);

  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);

  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition(int nbrParticleSector, int lzSector, int* nbrParticlesPerFlavorSector, RealVector& groundState, AbstractArchitecture* architecture);
};

#endif
