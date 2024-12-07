////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                            sign precalculation table                       //
//                                                                            //
//                        last modification : 12/12/2005                      //
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


#ifndef FERMIONONSPHEREWITHSPINALLSZ_H
#define FERMIONONSPHEREWITHSPINALLSZ_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereWithSpinAllSz :  public FermionOnSphereWithSpin
{

  friend class FermionOnSphereWithSpinAllSzLzSymmetry;
  friend class FermionOnSphereWithSpinAllSzLzSzSymmetry; 
  friend class FermionOnSphereWithSpinAllSzSzSymmetry;
  friend class FermionOnSphereWithSpinAllSzGutzwillerProjection;
/*   friend class FermionOnSphereWithSpinAllSzLzSymmetry; */

 protected:


  // temporary storage during state generation
  int **MaxTotalLz;


 public:

  // default constructor
  //
  FermionOnSphereWithSpinAllSz();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinAllSz (int nbrFermions, int totalLz, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinAllSz(const FermionOnSphereWithSpinAllSz& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinAllSz ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinAllSz& operator = (const FermionOnSphereWithSpinAllSz& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();


  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // Project the state from the tunneling space (all Sz's)
  // to the space with the fixed projection of Sz (given by SzValue)
  //
  // state = state that needs to be projected
  // su2Space = the subspace onto which the projection is carried out
  // SzValue = the desired value of Sz

  virtual RealVector ForgeSU2FromTunneling(RealVector& state, FermionOnSphereWithSpin& su2Space, int SzValue);

  // Project the state from the tunneling space (all Sz's)
  // to the U(1) space (u1Space)
  //
  // state = state that needs to be projected
  // u1Space = the subspace onto which the projection is carried out
  virtual RealVector ForgeU1FromTunneling(RealVector& state, FermionOnSphere& u1Space);

  // Calculate mean value <Sx> in a given state
  //
  // state = given state
  virtual double MeanSxValue(RealVector& state);

  // Calculate mean value <Sz> in a given state
  //
  // state = given state
  virtual double MeanSzValue(RealVector& state);

  // Artificially extend a state of a U(1) Hilbert space to a SU(2) space with all sz sectors
  //
  // state = state that needs to be projected
  // u1space = U(1) space of the input state
  // return value = input state expression in the SU(2) basis
  virtual ComplexVector U1ToSU2AllSz(ComplexVector& state, FermionOnSphere& u1space);

  // convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
  //
  // state = state that needs to be projected
  // su2space = SU(2) space with fixed sz of the input state
  // return value = input state expression in the SU(2) basis
  virtual RealVector SU2ToSU2AllSz(RealVector& state, ParticleOnSphereWithSpin* su2space);

  // Artificially extend a state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
  //
  // state = state that needs to be projected
  // su2space = SU(2) space with fixed sz of the input state
  // return value = input state expression in the SU(2) basis
  virtual ComplexVector SU2ToSU2AllSz(ComplexVector& state, FermionOnSphereWithSpin& su2space);

  // convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(RealVector& initialState, RealVector& targetState, RealMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

  // convert a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

  // recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  virtual void TransformOneBodyBasisRecursive(RealVector& targetState, double coefficient,
					      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, RealMatrix* oneBodyBasis);
  
  // recursive part of the convertion from a SU(2) basis to another one, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
				      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, ComplexMatrix* oneBodyBasis);


  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, 
									 bool removeBinomialCoefficient = false, AbstractArchitecture* architecture = 0);
   
  // core part of the entanglement matrix evaluation for the particle partition
  // 
  // minIndex = first index to consider in the complementary Hilbert space
  // nbrIndex = number of indices to consider in the complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // entanglementMatrix = pointer to entanglement matrix
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = number of components that have been added to the entanglement matrix
  virtual long EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, 
								       ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealMatrix* entanglementMatrix, 
								       bool removeBinomialCoefficient = false);
   
  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
  // weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
  // weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual RealMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
													 int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
													 int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, 
													 RealMatrix& entanglementMatrix);


 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);


  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up ( omitted: int totalSpin)
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, long pos);
  
};

#endif


