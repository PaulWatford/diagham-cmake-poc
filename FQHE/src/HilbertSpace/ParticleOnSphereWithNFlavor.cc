////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithNFlavor.h"

#include <iostream>
#include <cstring>

using namespace std;
using std::cout;
using std::endl;



/*****************************************************************/
/*  Destructor                                                   */
/*****************************************************************/

ParticleOnSphereWithNFlavor::~ParticleOnSphereWithNFlavor()
{
}


/*****************************************************************/
/*  Fermion statistic                                            */
/*****************************************************************/

/*int ParticleOnSphereWithNFlavor::GetParticleStatistic()
{
  return 1; // fermions
}
*/

/**************************************************************/
/*  Default ProdAdProdA (same behavior as SU4)               */
/**************************************************************/
// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
int ParticleOnSphereWithNFlavor::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

/*****************************************************************/
/*  Sum_sigma a^+_{m sigma} a_{m sigma}                         */
/*****************************************************************/

double ParticleOnSphereWithNFlavor::AdA(int index, int m)
{
  double result = 0.0;

  for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    result += this->AdsigmaAsigma(index, m, sigma);

  return result;
}

double ParticleOnSphereWithNFlavor::AdA(long index, int m)
{
  double result = 0.0;

  for (int sigma = 0; sigma < this->NbrFlavors; ++sigma)
    result += this->AdsigmaAsigma((int)index, m, sigma);

  return result;
}


/**************************************************************/
/*  Wave function wrappers (identical to SU4 pattern)        */
/**************************************************************/

Complex ParticleOnSphereWithNFlavor::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}


Complex ParticleOnSphereWithNFlavor::EvaluateWaveFunctionWithTimeCoherence (RealVector& state,
         RealVector& position, AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, this->HilbertSpaceDimension);
}

Complex ParticleOnSphereWithNFlavor::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
    int firstComponent, int nbrComponent)
{
return Complex(0.0, 0.0);
}

Complex ParticleOnSphereWithNFlavor::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
    AbstractFunctionBasis& basis, 
    int nextCoordinates, int firstComponent, 
    int nbrComponent)
{
return Complex(0.0, 0.0);
}

void ParticleOnSphereWithNFlavor::InitializeWaveFunctionEvaluation(bool timeCoherence)
{
  // nothing to do by default
}
/**************************************************************/
/*  Default density matrix (same as SU4 dummy)               */
/**************************************************************/

RealSymmetricMatrix ParticleOnSphereWithNFlavor::EvaluatePartialDensityMatrixParticlePartition(int nbrParticleSector, int lzSector, int* nbrParticlesPerFlavorSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}

