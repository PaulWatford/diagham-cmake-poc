#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "MathTools/BinomialCoefficients.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

// evaluate Hilbert space dimension for fermions within a single band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a two band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a two band and a cap on the number of particles per band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// maxBand0 = maximum number of particles in band 0
// maxBand1 = maximum number of particles in band 1
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentBand0 = current number of particles in band 0
// currentBand1 = current number of particles in band 1
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int maxBand0, int maxBand1, int currentKx, int currentKy, int currentBand0 = 0, int currentBand1 = 0, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a two band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a three band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionThreeBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a three band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a three band and spin and valley conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// nbrValleyPlus = number of particles with a valley plus  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionThreeBandWithSpinAndValleyEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int nbrValleyPlus, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a four band basis
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionFourBandEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a four band spin-conserved basis
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionFourBandWithSpinEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a four band spin-conserved basis and a cap on the number of particles per band
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// maxBand0 = maximum number of particles in band 0
// maxBand1 = maximum number of particles in band 1
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentBand0 = current number of particles in band 0
// currentBand1 = current number of particles in band 1
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int maxBand0, int maxBand1, int currentKx, int currentKy, int nbrSpinUp, int currentBand0 = 0, int currentBand1 = 0, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for bosions within a single band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long BosonSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for bosons within two band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long BosonTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for bosons within two band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long BosonTwoBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for bosons within three band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long BosonThreeBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions on a cubic lattice within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long FermionCubicLatticeSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);

// evaluate Hilbert space dimension for fermions on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);

// evaluate Hilbert space dimension for fermions on a cubic lattice within four bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);

// evaluate Hilbert space dimension for bosons on a cubic lattice within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long BosonCubicLatticeSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);

// evaluate Hilbert space dimension for bosons on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long BosonCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);

// evaluate Hilbert space dimension for bosons on a cubic lattice within four bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long BosonCubicLatticeFourBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);

// evaluate Hilbert space dimension for fermions on a hypercubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// nbrSiteT = number of sites along t
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// return value = Hilbert space dimension
long FermionHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0, int currentTotalKt = 0);

// evaluate Hilbert space dimension for bosons on a hypercubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// nbrSiteT = number of sites along t
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// return value = Hilbert space dimension
long BosonHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0, int currentTotalKt = 0);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETopInsulatorGetDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('t', "nbr-sitet", "number of sites along the t direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "bosons", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "no-inversion", "do not assume inversion symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "spin-conserved", "assume that the spin is conserved in the two band model");
  (*SystemGroup) += new BooleanOption  ('\n', "add-valley", "add a conserved  U(1) valley degree of freedom");
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption  ('\n', "4d", "consider a 4d model instead of a 2d model");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-band0", "maximum number of particles in band 0 (negative if this number should be equal to the number of orbitals)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-band1", "maximum number of particles in band 1 (negative if this number should be equal to the number of orbitals)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-band2", "maximum number of particles in band 2 (negative if this number should be equal to the number of orbitals)", -1);
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_topinsulator_nbrsubbands_n_nbrparticles_x_nbrsitex_y_nbrsitey.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorGetDimension -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  int NbrSitesZ = Manager.GetInteger("nbr-sitez"); 
  int NbrSitesT = Manager.GetInteger("nbr-sitet"); 
  
  if (Manager.GetInteger("nbr-subbands") == 1)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      if (Manager.GetBoolean("3d") == false)
		{
                  if ((Manager.GetBoolean("no-inversion") == true) || ((kx <= ((NbrSitesX - kx) % NbrSitesX)) && (ky <= ((NbrSitesY - ky) % NbrSitesY))))
                    {
                      long Dimension = 0;
                      if (Manager.GetBoolean("bosons") == false)
                        Dimension = FermionSingleBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
                      else
                        Dimension = BosonSingleBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
                      cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
                    }
                }
              else
                {
		  for (int kz = 0; kz < NbrSitesZ; ++kz)
		    {
                      if ((Manager.GetBoolean("no-inversion") == true) || ((kx <= ((NbrSitesX - kx) % NbrSitesX)) && (ky <= ((NbrSitesY - ky) % NbrSitesY)) && (kz <= ((NbrSitesZ - kz) % NbrSitesZ))))
                        {
                          long Dimension  = 0l;
                          if (Manager.GetBoolean("bosons") == false)
                              Dimension = FermionCubicLatticeSingleBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
                          else
                              Dimension = BosonCubicLatticeSingleBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
                          cout << "(kx=" << kx << ",ky=" << ky << ",kz=" << kz << ") : " << Dimension << endl;
                        }
		    }
                }
	    }
	}
    }
  if (Manager.GetInteger("nbr-subbands") == 2)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      if (Manager.GetBoolean("3d") == false)
		{
		  if (Manager.GetBoolean("4d") == false)
		    {
		      if ((Manager.GetBoolean("no-inversion") == true) || ((kx <= ((NbrSitesX - kx) % NbrSitesX)) && (ky <= ((NbrSitesY - ky) % NbrSitesY))))
			{
			  if (Manager.GetBoolean("spin-conserved") == false)
			    {
			      long Dimension = 0l;
			      if (Manager.GetBoolean("bosons") == false)
				{
				  if ((Manager.GetInteger("max-band0") < 0) && (Manager.GetInteger("max-band1") < 0))
				    {
				      Dimension = FermionTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
				    }
				  else
				    {
				      int MaxBand0 = Manager.GetInteger("max-band0");
				      int MaxBand1 = Manager.GetInteger("max-band1");
				      if (MaxBand0 < 0)
					{
					  MaxBand0 = NbrSitesX * NbrSitesY;
					}
				      if (MaxBand1 < 0)
					{
					  MaxBand1 = NbrSitesX * NbrSitesY;
					}
				      Dimension = FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, NbrSitesX - 1, NbrSitesY - 1);
				    }
				}
			      else
				{
				  Dimension = BosonTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
				}
			      cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
			    }
			  else
			    {
			      long TotalDimension = 0l;
			      for (int NbrSpinUp = 0; NbrSpinUp <= NbrParticles; ++NbrSpinUp)
				{
				  long Dimension = 0l;
				  if (Manager.GetBoolean("bosons") == false)
				    Dimension = FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp);
				  else
				    Dimension = BosonTwoBandWithSpinEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp);
				  TotalDimension += Dimension;
				  cout << "(kx=" << kx << ",ky=" << ky << ") 2Sz=" << ((2 * NbrSpinUp) - NbrParticles) << " : " << Dimension << endl;			  
				}
			      //		      cout << "(kx=" << kx << ",ky=" << ky << ") : " << TotalDimension << endl;
			    }
			}
		    }
		  else
		    {
		      long TotalDimension = 0l;
		      for (int kz = 0; kz < NbrSitesZ; ++kz)
			{
			  for (int kt = 0; kt < NbrSitesT; ++kt)
			    {
			      long Dimension = 0l;
			      if (Manager.GetBoolean("bosons") == false)
				Dimension = FermionHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, kt, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesT, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1, NbrSitesT - 1);
			      else
				Dimension = BosonHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, kt, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesT, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1, NbrSitesT - 1);
			      TotalDimension += Dimension;
			      cout << "(kx=" << kx << ",ky=" << ky << ",kz=" << kz << ",kt=" << kt << ") : " << Dimension << endl;
			    }
			}
		      cout << TotalDimension << endl;
		    }
		}
	      else
		{
		  long TotalDimension = 0l;
		  for (int kz = 0; kz < NbrSitesZ; ++kz)
		    {
		      long Dimension  = 0l;
		      if (Manager.GetBoolean("bosons") == false)
			Dimension = FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
		      else
			Dimension = BosonCubicLatticeTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
		      TotalDimension += Dimension;
		      cout << "(kx=" << kx << ",ky=" << ky << ",kz=" << kz << ") : " << Dimension << endl;
		    }
		  cout << TotalDimension << endl;
		}
	    }
	}
    }

  if (Manager.GetInteger("nbr-subbands") == 3)
    {
      long TotalDimension = 0l;
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      if (Manager.GetBoolean("3d") == false)
		{
		  if (Manager.GetBoolean("4d") == false)
		    {
		      if ((Manager.GetBoolean("no-inversion") == true) || 
			  ((kx <= ((NbrSitesX - kx) % NbrSitesX)) && (ky <= ((NbrSitesY - ky) % NbrSitesY))))
			{
			  if (Manager.GetBoolean("spin-conserved") == false)
			    {
			      long Dimension = 0l;
			      if (Manager.GetBoolean("bosons") == false)
				{
				  Dimension = FermionThreeBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
				}
			      else
				{
				  Dimension = BosonThreeBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
				}
			      TotalDimension += Dimension;
			      cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
			    }
			  else
			    {
			      if (Manager.GetBoolean("add-valley") == false)
				{
				  for (int NbrSpinUp = 0; NbrSpinUp <= NbrParticles; ++NbrSpinUp)
				    {
				      long Dimension = 0l;
				      if (Manager.GetBoolean("bosons") == false)
					{
					  Dimension = FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp);
					}
				      else
					{
					  cout << "warning : not implemented" << endl;
					  return 0;
					}
				      TotalDimension += Dimension;
				      cout << "(kx=" << kx << ",ky=" << ky << ") 2Sz=" << ((2 * NbrSpinUp) - NbrParticles) << " : " << Dimension << endl;			  
				    }
				}
			      else
				{
				  for (int NbrValleyPlus = 0; NbrValleyPlus <= NbrParticles; ++NbrValleyPlus)
				    {
				      for (int NbrSpinUp = 0; NbrSpinUp <= NbrParticles; ++NbrSpinUp)
					{
					  long Dimension = 0l;
					  if (Manager.GetBoolean("bosons") == false)
					    {
					      Dimension = FermionThreeBandWithSpinAndValleyEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp, NbrValleyPlus);
					    }
					  else
					    {
					      cout << "warning : not implemented" << endl;
					      return 0;
					    }
					  TotalDimension += Dimension;
					  cout << "(kx=" << kx << ",ky=" << ky << ") 2Pz=" << ((2 * NbrValleyPlus) - NbrParticles) << " 2Sz=" << ((2 * NbrSpinUp) - NbrParticles) << " : " << Dimension << endl;			  
					}
				    }
				}
			    }
			}
		    }
		  else
		    {
		      
		      cout << "warning : not implemented for 4d" << endl;
		      return 0;
		    }
		}
	      else
		{
		  cout << "warning : not implemented for 3d" << endl;
		  return 0;
		}
	    }
	}
      cout << TotalDimension << endl;
      // if (Manager.GetBoolean("spin-conserved") == false)
      // 	{
      // 	  BinomialCoefficients Binomial (3 * NbrSitesX * NbrSitesY);
      // 	  cout << Binomial(3 * NbrSitesX * NbrSitesY, NbrParticles) << endl;
      // 	}
      // else
      // 	{
      // 	  BinomialCoefficients Binomial (6 * NbrSitesX * NbrSitesY);
      // 	  cout << Binomial(6 * NbrSitesX * NbrSitesY, NbrParticles) << endl;
      // 	}
    }

  if (Manager.GetInteger("nbr-subbands") == 4)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      if (Manager.GetBoolean("3d") == true)
		{
		  long TotalDimension = 0l;
		  for (int kz = 0; kz < NbrSitesZ; ++kz)
		    {
		      long Dimension  = 0l;
		      if (Manager.GetBoolean("bosons") == false)
			Dimension = FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
		      else
			Dimension = BosonCubicLatticeFourBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
		      TotalDimension += Dimension;
		      cout << "(kx=" << kx << ",ky=" << ky << ",kz=" << kz << ") : " << Dimension << endl;
		    }
		  cout << TotalDimension << endl;
		}
	      else
		{
		  if (Manager.GetBoolean("spin-conserved") == false)
		    {
		      long Dimension  = 0l;
		      if (Manager.GetBoolean("bosons") == false)
			Dimension = FermionFourBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
		      cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
		    }
		  else
		    {
		      long TotalDimension = 0l;
		      for (int NbrSpinUp = 0; NbrSpinUp <= NbrParticles; ++NbrSpinUp)
			{
			  long Dimension = 0l;
			  if (Manager.GetBoolean("bosons") == false)
			    {
			      if ((Manager.GetInteger("max-band0") < 0) && (Manager.GetInteger("max-band1") < 0))
				{
				  Dimension = FermionFourBandWithSpinEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp);
				}
			      else
				{
				  int MaxBand0 = Manager.GetInteger("max-band0");
				  int MaxBand1 = Manager.GetInteger("max-band1");
				  if (MaxBand0 < 0)
				    {
				      MaxBand0 = 2 * NbrSitesX * NbrSitesY;
				    }
				  if (MaxBand1 < 0)
				    {
				      MaxBand1 = 2 * NbrSitesX * NbrSitesY;
				    }
				  Dimension = FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp);
				}
			    }
			  TotalDimension += Dimension;
			  cout << "(kx=" << kx << ",ky=" << ky << ") 2Sz=" << ((2 * NbrSpinUp) - NbrParticles) << " : " << Dimension << endl;
			}
		      cout << "(kx=" << kx << ",ky=" << ky << ") : " << TotalDimension << endl;
		    }
		}
	    }
	}
    }
      
  return 0;
}

// evaluate Hilbert space dimension for fermions within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// return value = Hilbert space dimension

long FermionSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    ++Count;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		++Count;
	    }
	}
      return Count;
    }
  Count += FermionSingleBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += FermionSingleBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a two band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// return value = Hilbert space dimension

long FermionTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  Count += FermionTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += (2 * FermionTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += FermionTwoBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}


// evaluate Hilbert space dimension for fermions within a two band and a cap on the number of particles per band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// maxBand0 = maximum number of particles in band 0
// maxBand1 = maximum number of particles in band 1
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentBand0 = current number of particles in band 0
// currentBand1 = current number of particles in band 1
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int maxBand0, int maxBand1, int currentKx, int currentKy, int currentBand0, int currentBand1, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((currentBand0 > maxBand0) || (currentBand1 > maxBand1))
    {
      return 0l;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    {
      return 0l;
    }
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  Count += FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentBand0 + 1, currentBand1 + 1);
  Count += FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentBand0 + 1, currentBand1);
  Count += FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentBand0, currentBand1 + 1);
  Count += FermionTwoBandWithBandCapEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentBand0, currentBand1);
  return Count;
}


// evaluate Hilbert space dimension for fermions within a two band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrParticles))
    return 0l;

  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count++;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count++;
	    }
	}
      return Count;
    }
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp, currentTotalKx, currentTotalKy);
  return Count;
}



// evaluate Hilbert space dimension for fermions within a three band basis
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionThreeBandEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 3l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 3l;
	    }
	}
      return Count;
    }
  Count += (FermionThreeBandEvaluateHilbertSpaceDimension(nbrFermions - 3,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (3 * FermionThreeBandEvaluateHilbertSpaceDimension(nbrFermions - 2,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (3 * FermionThreeBandEvaluateHilbertSpaceDimension(nbrFermions - 1,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += FermionThreeBandEvaluateHilbertSpaceDimension(nbrFermions,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a three band spin-conserved basis
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrFermionsUp, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermionsUp < 0) || (nbrFermionsUp > nbrFermions))
    return 0l;

  if (nbrFermions == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 3;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 3;
	    }
	}
      return Count;
    }
  Count += FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 6, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 3, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy));
  
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 5, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 3, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy));
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 5, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy));
  
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 4, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 3, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  Count += 9 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 4, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 4, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));


  Count += FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 3, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy));
  Count += 9 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy));
  Count += 9 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy));
  Count += FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy));

  
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += 9 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));

  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += 3 * FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx + currentKx, currentTotalKy + currentKy);

  Count += FermionThreeBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx, currentTotalKy);
  return Count;
}


// evaluate Hilbert space dimension for fermions within a three band and spin and valley conserved basis
//
// nbrFermions = number of nbrFermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of fermions with a spin up  
// nbrValleyPlus = number of fermions with a valley plus  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionThreeBandWithSpinAndValleyEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int nbrValleyPlus, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrFermions) || (nbrValleyPlus < 0) || (nbrValleyPlus > nbrFermions))
    return 0l;

  if (nbrFermions == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 3;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 3;
	    }
	}
      return Count;
    }
  for (int i = 12; i >= 0; --i)
    {
      int MaxS = 6;
       if (i < MaxS)
	{
	  MaxS = i;
	}
      for (int s = MaxS; s >= 0; --s)
	{	  
	  for (int k = MaxS; k >= 0; --k)
	    {
	      Count += FermionThreeBandWithSpinAndValleyEvaluateHilbertSpaceDimension(nbrFermions - i, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp - s, nbrValleyPlus - k, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
	    }
	}
    }
  return Count;
}

// evaluate Hilbert space dimension for fermions within a four band basis
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionFourBandEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 4l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 4l;
	    }
	}
      return Count;
    }
  Count += FermionFourBandEvaluateHilbertSpaceDimension(nbrFermions - 4,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  Count += (4 * FermionFourBandEvaluateHilbertSpaceDimension(nbrFermions - 3,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (6 * FermionFourBandEvaluateHilbertSpaceDimension(nbrFermions - 2,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (4 * FermionFourBandEvaluateHilbertSpaceDimension(nbrFermions - 1,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += FermionFourBandEvaluateHilbertSpaceDimension(nbrFermions,  kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a four band spin-conserved basis
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionFourBandWithSpinEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrFermionsUp, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)|| (nbrFermionsUp < 0) ||  (nbrFermionsUp > nbrFermions))
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  
  Count += FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 4, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  
  Count += (2 * FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (2 * FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  
  Count += FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 2, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += (2 * FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (2 * FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  
  Count += (2 * FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += (2 * FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  
  Count += FermionFourBandWithSpinEvaluateHilbertSpaceDimension(nbrFermions, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrFermionsUp, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a four band spin-conserved basis and a cap on the number of particles per band
//
// nbrFermions = number of fermions
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// maxBand0 = maximum number of particles in band 0
// maxBand1 = maximum number of particles in band 1
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentBand0 = current number of particles in band 0
// currentBand1 = current number of particles in band 1
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(int nbrFermions, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int maxBand0, int maxBand1, int currentKx, int currentKy, int nbrSpinUp, int currentBand0, int currentBand1, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)|| (nbrSpinUp < 0) ||  (nbrSpinUp > nbrFermions) || (currentBand0 > maxBand0) || (currentBand1 > maxBand1))
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      int TmpIncrement = 0;
      if (currentBand0 < maxBand0)
	{
	  TmpIncrement++;
	}
      if (currentBand1 < maxBand1)
	{
	  TmpIncrement++;
	}
      
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    {
	      Count += TmpIncrement;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		{
		  Count += TmpIncrement;
		}
	    }
	}
      return Count;
    }
  
  Count += FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 4, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 2, currentBand0 + 2, currentBand1 + 2, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 2, currentBand0 + 2, currentBand1 + 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 2, currentBand0 + 1, currentBand1 + 2, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0 + 2, currentBand1 + 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += ( FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 3, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0 + 1, currentBand1 + 2, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  
  Count += FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 2, currentBand0 + 1, currentBand1 + 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0 + 2, currentBand1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (2 * FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0 + 1, currentBand1 + 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0, currentBand1 + 2, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp, currentBand0 + 1, currentBand1 + 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp, currentBand0 + 1, currentBand1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp, currentBand0, currentBand1 + 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0 + 1, currentBand1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += (FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp - 1, currentBand0, currentBand1 + 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  
  Count += FermionFourBandWithSpinAndBandCapEvaluateHilbertSpaceDimension(nbrFermions, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, maxBand0, maxBand1, currentKx, currentKy - 1, nbrSpinUp, currentBand0, currentBand1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// return value = Hilbert space dimension

long BosonSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    ++Count;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		++Count;
	    }
	}
      return Count;
    }
  for (int i = nbrParticles; i >= 0; --i)
    Count += BosonSingleBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}

// evaluate Hilbert space dimension for bosons within two bands
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  for (int i = nbrParticles; i >= 0; --i)
    Count += ((long) i + 1l) * BosonTwoBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}

// evaluate Hilbert space dimension for bosons within two band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// nbrSpinUp = number of particles with a spin up  
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonTwoBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }

  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrParticles))
    return 0l;

  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count++;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count++;
	    }
	}
      return Count;
    }
  for (int i = nbrParticles; i >= 0; --i)
    for (int j = i; j >= 0; --j)
      Count += BosonTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp - j, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}

// evaluate Hilbert space dimension for bosons within three bands
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonThreeBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 3l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 3l;
	    }
	}
      return Count;
    }
  for (int i = nbrParticles; i >= 0; --i)
    Count += ((((long) i + 1l) * ((long) i + 2l)) / 2l) * BosonThreeBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}

// evaluate Hilbert space dimension for fermions on a cubic lattice within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long FermionCubicLatticeSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ,
        int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
    if (currentKz < 0)
    {
        currentKz = nbrSiteZ - 1;
        currentKy--;
        if (currentKy < 0)
        {
            currentKy = nbrSiteY - 1;
            currentKx--;
        }
    }
    if (nbrParticles == 0)
    {
        if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
                && ((currentTotalKz % nbrSiteZ) == kzMomentum))
            return 1l;
        else
            return 0l;
    }
    if (currentKx < 0)
        return 0l;
    long Count = 0;
    if (nbrParticles == 1)
    {
        for (int k = currentKz; k >= 0; --k)
        {
            if (((currentKx + currentTotalKx) % nbrSiteX == kxMomentum) && ((currentKy + currentTotalKy) % nbrSiteY == kyMomentum)
                    && ((k + currentTotalKz) % nbrSiteZ == kzMomentum))
                ++Count;
        }
        for (int j = currentKy - 1; j >= 0; --j)
        {
            for (int k = nbrSiteZ - 1; k >= 0; --k)
            {
                if (((currentKx + currentTotalKx) % nbrSiteX == kxMomentum) && ((j + currentTotalKy) % nbrSiteY == kyMomentum)
                        && ((k + currentTotalKz) % nbrSiteZ == kzMomentum))
                    ++Count;
            }
        }
        for (int i = currentKx - 1; i >=0; --i)
        {
            for (int j = nbrSiteY - 1; j >= 0; --j)
            {
                for (int k = nbrSiteZ - 1; k >= 0; --k)
                {
                    if (((i + currentTotalKx) % nbrSiteX == kxMomentum) && ((j + currentTotalKy) % nbrSiteY == kyMomentum)
                            && ((k + currentTotalKz) % nbrSiteZ == kzMomentum))
                        ++Count;
                }
            }
        }
        return Count;
    }
    Count += FermionCubicLatticeSingleBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ,
            currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz);
    Count += FermionCubicLatticeSingleBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ,
            currentKx, currentKy, currentKz - 1, currentTotalKx, currentTotalKy, currentTotalKz);
    return Count;
}

// evaluate Hilbert space dimension for fermions on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
  if (currentKz < 0)
    {
      currentKz = nbrSiteZ - 1;
      currentKy--;
      if (currentKy < 0)
	{
	  currentKy = nbrSiteY - 1;
	  currentKx--;
	}
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int k = currentKz; k >= 0; --k)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((currentKy + currentTotalKy) % nbrSiteY) == kyMomentum) && 
	      (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
	    Count += 2l;
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = nbrSiteZ - 1; k >= 0; --k)
	    {
	      if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		  && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		Count += 2l;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = nbrSiteZ - 1; k >= 0; --k)
		{
		  if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		      && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		    Count += 2l;
		}
	    }
	}
      return Count;
    }
  Count += FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz));
  Count += (2 * FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz));
  Count += FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx, currentTotalKy, currentTotalKz);
  return Count;
}

// evaluate Hilbert space dimension for bosons on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long BosonCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
  if (currentKz < 0)
    {
      currentKz = nbrSiteZ - 1;
      currentKy--;
      if (currentKy < 0)
	{
	  currentKy = nbrSiteY - 1;
	  currentKx--;
	}
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  for (int i = nbrParticles; i >= 0; --i)
    Count += (((long) i) + 1l) * BosonCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz));
  return Count;
}


// evaluate Hilbert space dimension for fermions on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
  if (currentKz < 0)
    {
      currentKz = nbrSiteZ - 1;
      currentKy--;
      if (currentKy < 0)
	{
	  currentKy = nbrSiteY - 1;
	  currentKx--;
	}
    }
  if (nbrParticles < 0)
    {
      return 0;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int k = currentKz; k >= 0; --k)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((currentKy + currentTotalKy) % nbrSiteY) == kyMomentum) && 
	      (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
	    Count += 4l;
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = nbrSiteZ - 1; k >= 0; --k)
	    {
	      if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		  && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		Count += 4l;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = nbrSiteZ - 1; k >= 0; --k)
		{
		  if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		      && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		    Count += 4l;
		}
	    }
	}
      return Count;
    }
  Count += FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(nbrParticles - 4, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), currentTotalKz + (4 * currentKz));
  Count += (4 * FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(nbrParticles - 3, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), currentTotalKz + (3 * currentKz)));
  Count += (6 * FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz)));
  Count += (4 * FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz));
  Count += FermionCubicLatticeFourBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx, currentTotalKy, currentTotalKz);
  return Count;
}

// evaluate Hilbert space dimension for bosons on a cubic lattice within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long BosonCubicLatticeSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
    if (currentKz < 0)
    {
        currentKz = nbrSiteZ - 1;
        currentKy--;
        if (currentKy < 0)
        {
            currentKy = nbrSiteY - 1;
            currentKx--;
        }
    }
    if (nbrParticles == 0)
    {
        if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
                && ((currentTotalKz % nbrSiteZ) == kzMomentum))
            return 1l;
        else
            return 0l;
    }
    if (currentKx < 0)
        return 0l;
    long Count = 0;
    for (int i = nbrParticles; i >= 0; --i)
        Count += BosonCubicLatticeSingleBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, kzMomentum,
                nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1,
                currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz));
    return Count;
}

// evaluate Hilbert space dimension for bosons on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long BosonCubicLatticeFourBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
  if (currentKz < 0)
    {
      currentKz = nbrSiteZ - 1;
      currentKy--;
      if (currentKy < 0)
	{
	  currentKy = nbrSiteY - 1;
	  currentKx--;
	}
    }
  if (nbrParticles < 0)
    {
      return 0;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int k = currentKz; k >= 0; --k)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((currentKy + currentTotalKy) % nbrSiteY) == kyMomentum) && 
	      (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
	    Count += 4l;
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = nbrSiteZ - 1; k >= 0; --k)
	    {
	      if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		  && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		Count += 4l;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = nbrSiteZ - 1; k >= 0; --k)
		{
		  if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		      && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		    Count += 4l;
		}
	    }
	}
      return Count;
    }
  for (int i = nbrParticles; i >= 0; --i)
    Count += ((((long) i + 1l) * ((long) i + 2l) * ((long) i + 3l)) / 6l) * BosonCubicLatticeFourBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz));
  return Count;
}

// evaluate Hilbert space dimension for fermions on a hypercubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// nbrSiteT = number of sites along t
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// return value = Hilbert space dimension

long FermionHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt)
{
  cout << nbrParticles << " " <<  kxMomentum << " " <<  kyMomentum << " " <<  kzMomentum << " " <<  ktMomentum << " " <<  nbrSiteX << " " <<  nbrSiteY << " " <<  nbrSiteZ << " " <<  nbrSiteT << " " <<  currentKx << " " <<  currentKy << " " <<  currentKz << " " <<  currentKt << " " <<  currentTotalKx << " " <<  currentTotalKy << " " <<  currentTotalKz << " " <<  currentTotalKt << endl;
  if (currentKt < 0)
    {
      currentKt = nbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = nbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = nbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum) && ((currentTotalKt % nbrSiteT) == ktMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int l = currentKt; l >= 0; --l)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((currentKy + currentTotalKy) % nbrSiteY) == kyMomentum) && 
	      (((currentKz + currentTotalKz) % nbrSiteZ) == kzMomentum) && (((l + currentTotalKt) % nbrSiteT) == ktMomentum))
	    Count += 2l;
	}
      for (int k = currentKz; k >= 0; --k)
	{
	  for (int l = nbrSiteT - 1; l >= 0; --l)
	    {
	      if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((currentKy + currentTotalKy) % nbrSiteY) == kyMomentum) && 
		  (((k + currentTotalKz) % nbrSiteZ) == kzMomentum) && (((l + currentTotalKt) % nbrSiteT) == ktMomentum))
		Count += 2l;
	    }
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = nbrSiteZ - 1; k >= 0; --k)
	    {
	      for (int l = nbrSiteT - 1; l >= 0; --l)
		{
		  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		      && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum) && (((l + currentTotalKt) % nbrSiteT) == ktMomentum))
		    Count += 2l;
		}
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = nbrSiteZ - 1; k >= 0; --k)
		{
		  for (int l = nbrSiteT - 1; l >= 0; --l)
		    {
		      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
			  && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum) && (((l + currentTotalKt) % nbrSiteT) == ktMomentum))
			Count += 2l;
		    }
		}
	    }
	}
      return Count;
    }
  Count += FermionHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, kzMomentum, ktMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, nbrSiteT, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), currentTotalKt + (2 * currentKt));
  Count += (2 * FermionHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, kzMomentum, ktMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, nbrSiteT, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt));
  Count += FermionHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, kzMomentum, ktMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, nbrSiteT, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx, currentTotalKy, currentTotalKz, currentTotalKt);
  return Count;
}


// evaluate Hilbert space dimension for bosons on a hypercubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// nbrSiteT = number of sites along t
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// return value = Hilbert space dimension

long BosonHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt)
{
  if (currentKt < 0)
    {
      currentKt = nbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = nbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = nbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum)&& ((currentTotalKt % nbrSiteT) == ktMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  for (int i = nbrParticles; i >= 0; --i)
    Count += (((long) i) + 1l) * BosonHyperCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - i, kxMomentum, kyMomentum, kzMomentum, ktMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, nbrSiteT, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz), currentTotalKt + (i * currentKt));
  return Count;
}
