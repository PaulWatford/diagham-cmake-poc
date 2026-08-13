#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithPolarizedSpin.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
// #include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry.h"
// #include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry.h"
// #include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpace.h"
// #include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareOpenLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeFromFileInteractionOneBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian.h"
 
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBandStructure.h"
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBlochHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FTIGenericInteractionFromFileSingleBand" , "0.02");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "open-y", "only enforce momentum conservation x, treating ky as a coordinate along y");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "name of the file containing the two-body interaction matrix elements");
  (*SystemGroup) += new SingleStringOption  ('\n', "threebody-interaction-file", "name of the file containing an optional three-body interaction matrix elements");
  //  (*SystemGroup) += new BooleanOption  ('\n', "real-interaction", "assume that the two-body interaction matrix elements are real");
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-name", "name of the two-body interaction", "noname");
  (*SystemGroup) += new SingleDoubleOption  ('u', "interaction-rescaling", "global rescaling of the two-body interactio", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('u', "threebody-interaction-rescaling", "global rescaling of the three-body interaction", 1.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "singleparticle-file", "optional name of the file containing the one-body matrix elements");
  //  (*SystemGroup) += new BooleanOption  ('\n', "full-singleparticle", "the one-body matrix element file contains off-diagonal inter-band contributions");
  //  (*SystemGroup) += new BooleanOption  ('\n', "complex-singlebody", "the one-body matrix element file contains complex entries (only valid when using --full-singleparticle)");
  (*SystemGroup) += new BooleanOption  ('\n', "add-valley", "add valley-like degree of freedom (i.e. U(1) symmetry) included in --interaction-file");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pz-value", "twice the valley Pz value", 0);
  //  (*SystemGroup) += new SingleIntegerOption  ('\n', "ez-value", "twice the Ez =1/2(N_{1u}+N_{2d}-N_{1d}-N_{2u}) value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "add-spin", "add spin 1/2 degree of freedom while assuming an SU(2) invariant interaction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-value", "twice the spin Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-valleyspin", "use the spin per valley instead of --sz-value and --ez-value");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz1-value", "twice the Sz value in valley 1", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz2-value", "twice the Sz value in valley 2", 0);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-sectors", "provide an ascii file that indicates which symmetry sectors have to be computed");
  (*SystemGroup) += new BooleanOption ('\n', "variable-nbreigenvalues", "when using the --selected-points option, use the last column of the momentum file to fix the number of eigenvalues to compute (overriding the command line options)", false);
  (*SystemGroup) += new BooleanOption  ('\n', "disable-pzsymmetry", "disable the valley Pz<->-Pz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "disable-szsymmetry", "disable the valley Sz<->-Sz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "apply a temporary energy shift during the diagonalization", 0.0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption ('\n', "hilbert-directory", "directory where Hilbert spaces should be read and stored (bypassing their construction)");
  (*PrecalculationGroup) += new BooleanOption ('\n', "hilbert-only", "only compute the Hilbert spaces, bypassing the diagonalization");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hilbert", "show the Hilbert space basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIGenericInteractionFromFileSingleBand -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("interaction-file") == 0) && (Manager.GetString("threebody-interaction-file") == 0))
    {
      cout << "no interaction file defined" << endl;
      cout << "see man page for option syntax or type FTIGenericInteractionFromFileSingleBand -h" << endl;
      return -1;
    }
  if ((Manager.GetString("interaction-file") != 0) && (!(IsFile(Manager.GetString("interaction-file")))))
    {
      cout << "interaction file " << Manager.GetString("interaction-file")<< " does not exist" << endl;
      return -1;
    }
  if ((Manager.GetString("threebody-interaction-file") != 0) && (!(IsFile(Manager.GetString("threebody-interaction-file")))))
    {
      cout << "three-body interaction file " << Manager.GetString("threebody-interaction-file")<< " does not exist" << endl;
      return -1;
    }

  //  if ((Manager.GetBoolean("real-interaction") == true) && (Manager.GetBoolean("conserve-bandoccuption") == true))
  //    {
  //      cout << "warning, be sure that your interaction preserves band occupation when using --real-interaction" << endl;
  //    }
    
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSites = 2 * NbrSitesX * NbrSitesY;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  double EnergyShift = Manager.GetDouble("energy-shift");
  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  if (Manager.GetBoolean("open-y") == true)
    {
      MinKy = 0;
      MaxKy = 0;
    }
  int NbrMomentumSectors = (MaxKx - MinKx + 1) * (MaxKy - MinKy + 1);
  int MinSz = NbrParticles & 1;
  int MaxSz = MinSz;
  if (Manager.GetBoolean("add-spin") == true)
    {
      MinSz = Manager.GetInteger("sz-value") | (NbrParticles & 1);
      if (Manager.GetBoolean("use-valleyspin") == true)
	{
	  MinSz = (Manager.GetInteger("sz1-value") + Manager.GetInteger("sz2-value")) | (NbrParticles & 1);
	}
      MaxSz = MinSz;
    }  
  int MinPz = NbrParticles & 1;
  int MaxPz = MinPz;
  if (Manager.GetBoolean("add-valley") == true)
    {
      MinPz = Manager.GetInteger("pz-value") | (NbrParticles & 1);
      MaxPz = MinPz;
    }  
  // int MinEz = NbrParticles & 1;
  // int MaxEz = MinEz;
  // if ((Manager.GetBoolean("add-valley") == true) && (Manager.GetBoolean("add-spin") == true))
  //   {
  //     MinEz = Manager.GetInteger("ez-value") | (NbrParticles & 1);
  //     if (Manager.GetBoolean("use-valleyspin") == true)
  // 	{
  // 	  MinEz = (Manager.GetInteger("sz1-value") - Manager.GetInteger("sz2-value")) | (NbrParticles & 1);
  // 	}
  //     MaxEz = MinEz;
  //   }  
  bool DisablePzMinusPzSymmetry = Manager.GetBoolean("disable-pzsymmetry");
  bool DisableSzMinusSzSymmetry = Manager.GetBoolean("disable-szsymmetry");
  bool UsePzMinusPzSymmetry = !DisablePzMinusPzSymmetry;
  bool UseSzMinusSzSymmetry = !DisableSzMinusSzSymmetry;
    
  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }


  char* FileSystemGeometry = new char [512];
  char* CommentLine = new char [256];
  if (Manager.GetBoolean("add-valley") == false)
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  //	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	      if (Manager.GetBoolean("open-y") == false)
		{
		  sprintf (CommentLine, "eigenvalues\n# kx ky");
		}
	      else
		{
		  sprintf (CommentLine, "eigenvalues\n# kx");
		}
	    }
	  // else
	  //   {
	  //     sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	  //     sprintf (CommentLine, "eigenvalues\n# kx ky n1 n2");
	  //   }
	}
      else
	{
	  // if (Manager.GetBoolean("conserve-bandoccuption") == false)
	     {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	      if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		{
		  UseSzMinusSzSymmetry = true;
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (CommentLine, "eigenvalues\n# Sz Szsym kx ky");
		    }
		  else
		    {
		      sprintf (CommentLine, "eigenvalues\n# Sz Szsym kx");
		    }		  
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (CommentLine, "eigenvalues\n# Sz kx ky");
		    }
		  else
		    {
		      sprintf (CommentLine, "eigenvalues\n# Sz kx");
		    }
		}
	    }
	  // else
	  //   {
	  //     sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	  //     sprintf (CommentLine, "eigenvalues\n# Sz kx ky n1up n1down n2up n2down");
	  //   }
	}
   }
  else
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  //	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz);
	      if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		{
		  UsePzMinusPzSymmetry = true;
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz Pzsym kx ky");
		    }
		  else
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz Pzsym kx");
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz kx ky");
		    }
		  else
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz kx");
		    }
		}
	    }
	  // else
	  //   {
	  //     sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz);
	  //     sprintf (CommentLine, "eigenvalues\n# Pz kx ky n1plus n1minus n2plus n2minus");
	  //   }
	}
      else
	{
	  // if (Manager.GetBoolean("conserve-bandoccuption") == false)
	  {
	       sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz, MinSz);
	      if (((MinSz == 0) && (DisableSzMinusSzSymmetry == false)) || ((MinPz == 0) && (DisablePzMinusPzSymmetry == false)))
		{
		  UseSzMinusSzSymmetry = true;
		  UsePzMinusPzSymmetry = true;
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz Sz Pzsym Szsym kx ky");
		    }
		  else
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz Sz Pzsym Szsym kx");
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz Sz kx ky");
		    }
		  else
		    {
		      sprintf (CommentLine, "eigenvalues\n# Pz Sz kx");
		    }
		}
	    }
	  // else
	  //   {
	  //     sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d_ez_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz, MinEz, MinSz);
	  //     sprintf (CommentLine, "eigenvalues\n# Pz Sz Ez kx ky n1upplus n2upplus n1upminus n2upminus n1downplus n2downplus n1downminus n2downminus");
	  //   }
	}
    }
  char* BandCapPrefix = 0;
  BandCapPrefix = new char [2];
  sprintf (BandCapPrefix, "");

  char* FilePrefix = new char [512 + strlen(FileSystemGeometry) + strlen(BandCapPrefix)+ strlen(Manager.GetString("interaction-name"))];
  if (Manager.GetBoolean("flat-band"))
    {
      if (Manager.GetBoolean("open-y") == false)
	{
	  sprintf (FilePrefix, "%s_oneband_flatband_%s%s_%s", StatisticPrefix, Manager.GetString("interaction-name"), BandCapPrefix, FileSystemGeometry);
	}
      else
	{
	  sprintf (FilePrefix, "%s_oneband_openy_flatband_%s%s_%s", StatisticPrefix, Manager.GetString("interaction-name"), BandCapPrefix, FileSystemGeometry);
	}
    }
  else
    {
      if (Manager.GetBoolean("open-y") == false)
	{
	  sprintf (FilePrefix, "%s_oneband_u_%.3f_%s%s_%s", StatisticPrefix, Manager.GetDouble("interaction-rescaling"),
		   Manager.GetString("interaction-name"), BandCapPrefix, FileSystemGeometry);
	}
      else
	{
	  sprintf (FilePrefix, "%s_oneband_openy_u_%.3f_%s%s_%s", StatisticPrefix, Manager.GetDouble("interaction-rescaling"),
		   Manager.GetString("interaction-name"), BandCapPrefix, FileSystemGeometry);
	}
    }
  
  char* EigenvalueOutputFile = new char [512 + strlen(FilePrefix)];
  
  if (Manager.GetString("eigenvalue-file") != 0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    }
  else
    {
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);
    }
  

  Abstract2DTightBindingModel* TightBindingModel;
  
  if (Manager.GetString("import-onebody") == 0)
    {
      if (Manager.GetString("singleparticle-file") == 0)
	{
	  double* DummyChemicalPotentials = new double[1];
	  DummyChemicalPotentials[0] = 0.0;
	  TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, 1, DummyChemicalPotentials,
									 0.0, 0.0, Architecture.GetArchitecture(), true);
	}
      else
	{
	  double* DummyChemicalPotentials = new double[1];
	  DummyChemicalPotentials[0] = 0.0;
	  int* TmpKxValues = 0;
	  int* TmpKyValues = 0;
	  int* TmpValleyIndices = 0;
	  int* TmpBandIndices1 = 0;	  
	  int* TmpBandIndices2 = 0;	  
	  double* TmpOneBodyEnergies = 0;
	  MultiColumnASCIIFile OneBodyEnergyFile;
	  if (OneBodyEnergyFile.Parse(Manager.GetString("singleparticle-file")) == false)
	    {
	      OneBodyEnergyFile.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if (OneBodyEnergyFile.GetNbrLines() == 0)
	    {
	      cout << Manager.GetString("singleparticle-file") << " is an empty file" << endl;
	      return 0;
	    }
	  int NbrEnergies = OneBodyEnergyFile.GetNbrLines();
	  if (((Manager.GetBoolean("add-valley") == false) && (NbrEnergies != (NbrSitesX * NbrSitesY)))
	      || ((Manager.GetBoolean("add-valley") == true) && (NbrEnergies != (2 * NbrSitesX * NbrSitesY))))
	    {
	      cout << Manager.GetString("singleparticle-file") << " has a wrong number of lines (has "
		   << NbrEnergies << ", should be " << (NbrSitesX * NbrSitesY)
		   << " without valley, " << (2 * NbrSitesX * NbrSitesY) << " with valley)" << endl;
	      return 0;
	    }
	  if (((Manager.GetBoolean("add-valley") == false) && (OneBodyEnergyFile.GetNbrColumns() < 3))
	      || ((Manager.GetBoolean("add-valley") == true) && (OneBodyEnergyFile.GetNbrColumns() < 4)))
	    {
	      cout << Manager.GetString("singleparticle-file") << " has a wrong number of columns (has "
		   << OneBodyEnergyFile.GetNbrColumns() << ", should be at least 3 without valley, 4 with valley)" << endl;
	      return 0;
	    }

	  int TmpNbrBands = 1;
	  if (Manager.GetBoolean("add-valley") == true)
	    {
	      TmpNbrBands = 2;
	    }
	  if (Manager.GetBoolean("add-spin") == true)
	    {
	      TmpNbrBands *= 2;
	    }
	  int* KxValues = new int[NbrSitesX * NbrSitesY];
	  int* KyValues = new int[NbrSitesX * NbrSitesY];
	  int TmpIndex = 0;
	  int* IndexToKIndex = new int[NbrSitesX * NbrSitesY];
	  int* NbrBandsPerKSector = new int[NbrSitesX * NbrSitesY];
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
		{
		  KxValues[TmpIndex] = i;
		  KyValues[TmpIndex] = j;
		  IndexToKIndex[(i * NbrSitesY) + j] = TmpIndex;
		  ++TmpIndex;
		}
	    }
	  TmpKxValues = OneBodyEnergyFile.GetAsIntegerArray(0);
	  TmpKyValues = OneBodyEnergyFile.GetAsIntegerArray(1);

	  double** Energies = new double*[NbrSitesX * NbrSitesY];
	  TmpIndex = 0;
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
		{
		  Energies[TmpIndex] = new double[TmpNbrBands];
		  ++TmpIndex;
		}
	    }
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(2);
	      if (TmpOneBodyEnergies == 0)
		{
		  OneBodyEnergyFile.DumpErrors(cout) << endl;
		  return 0;
		}
	      
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  for (int i = 0; i < NbrEnergies; ++i)
		    {
		      Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][0] = TmpOneBodyEnergies[i];
		    }
		}
	      else
		{
		  for (int i = 0; i < NbrEnergies; ++i)
		    {
		      Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][0] = TmpOneBodyEnergies[i];
		      Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][1] = TmpOneBodyEnergies[i];
		    }
		}			      
	    }
	  else
	    {
	      TmpValleyIndices = OneBodyEnergyFile.GetAsIntegerArray(2);
	      TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(3);
	      if (TmpOneBodyEnergies == 0)
		{
		  OneBodyEnergyFile.DumpErrors(cout) << endl;
		  return 0;
		}
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  for (int i = 0; i < NbrEnergies; ++i)
		    {
		      Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][((TmpValleyIndices[i] + 1) >> 1)] = TmpOneBodyEnergies[i];
			}
		}
	      else
		{
		  for (int i = 0; i < NbrEnergies; ++i)
		    {
		      Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][(((TmpValleyIndices[i] + 1) >> 1))] = TmpOneBodyEnergies[i];
		      Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][2 + ((TmpValleyIndices[i] + 1) >> 1)] = TmpOneBodyEnergies[i];
		    }
		}			      
	    }
	  TightBindingModel = new TightBindingModel2DExplicitBandStructure (NbrSitesX, NbrSitesY, TmpNbrBands, 0.0, 0.0,
									    KxValues, KyValues, Energies,
									    Architecture.GetArchitecture(), true);
	  TmpIndex = 0;
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
		{
		  delete[] Energies[TmpIndex];
		  ++TmpIndex;
		}
	    }
	  delete[] Energies;
	  char* BandStructureOutputFile = new char [64 + strlen(FilePrefix)];
	  sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  TightBindingModel->WriteBandStructure(BandStructureOutputFile);

	  delete[] IndexToKIndex;
	  delete[] KxValues;
	  delete[] KyValues;
	}      
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }

  bool FirstRunFlag = true;

  // if (Manager.GetBoolean("real-interaction"))
  //   {
  //     Lanczos.SetRealAlgorithms();
  //   }

  int* NbrParticlesBand1UpPlus = 0;
  int* NbrParticlesBand1UpMinus = 0;
  int* NbrParticlesBand1DownPlus = 0;
  int* NbrParticlesBand1DownMinus = 0;
  int* KxMomenta = 0;
  int* KyMomenta = 0;
  int* SzValues = 0;
  int* PzValues = 0;
  int* EzValues = 0;
  int* SzParityValues1 = 0;
  int* SzParityValues2 = 0;
  int* PzParityValues1 = 0;
  int* PzParityValues2 = 0;
  int NbrSymmetrySectors = NbrMomentumSectors;
  int* NbrRequestedEigenstates = 0;
  if (Manager.GetString("selected-sectors") == 0)
    {
      //      if (Manager.GetBoolean("conserve-bandoccuption") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  // no spin, no valley
		  NbrSymmetrySectors = 1;
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = NbrParticles;
		  NbrParticlesBand1UpMinus[0] = 0;
		  NbrParticlesBand1DownPlus[0] = 0;
		  NbrParticlesBand1DownMinus[0] = 0;
		  SzParityValues1[0] = 0;
		  SzParityValues2[0] = 0;
		  PzParityValues1[0] = 0;
		  PzParityValues2[0] = 0;
		}
	      else
		{
		  // valley, no spin
		  NbrSymmetrySectors = 1;
		  if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		    {
		      NbrSymmetrySectors = 2;
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = (NbrParticles + MinPz) / 2;
		  NbrParticlesBand1UpMinus[0] = (NbrParticles - MinPz) / 2;
		  NbrParticlesBand1DownPlus[0] = 0;
		  NbrParticlesBand1DownMinus[0] = 0;
		  SzParityValues1[0] = 0;
		  SzParityValues2[0] = 0;
		  PzParityValues2[0] = 0;
		  if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		    {
		      PzParityValues1[0] = 1;
		      NbrParticlesBand1UpPlus[NbrMomentumSectors] = (NbrParticles + MinPz) / 2;
		      NbrParticlesBand1UpMinus[NbrMomentumSectors] = (NbrParticles - MinPz) / 2;
		      NbrParticlesBand1DownPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1DownMinus[NbrMomentumSectors] = 0;
		      SzParityValues1[NbrMomentumSectors] = 0;
		      PzParityValues1[NbrMomentumSectors] = -1;
		      SzParityValues2[NbrMomentumSectors] = 0;
		      PzParityValues2[NbrMomentumSectors] = 0;
		    }
		  else
		    {
		      PzParityValues1[0] = 0;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  // spin, no valley
		  NbrSymmetrySectors = 1;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      NbrSymmetrySectors = 2;
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz) / 2;
		  NbrParticlesBand1UpMinus[0] = 0;
		  NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz) / 2;
		  NbrParticlesBand1DownMinus[0] = 0;
		  PzParityValues1[0] = 0;
		  PzParityValues2[0] = 0;
		  SzParityValues2[0] = 0;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      SzParityValues1[0] = 1;
		      NbrParticlesBand1UpPlus[NbrMomentumSectors] = (NbrParticles + MinPz) / 2;
		      NbrParticlesBand1UpMinus[NbrMomentumSectors] = (NbrParticles - MinPz) / 2;
		      NbrParticlesBand1DownPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1DownMinus[NbrMomentumSectors] = 0;
		      PzParityValues1[NbrMomentumSectors] = 0;
		      SzParityValues1[NbrMomentumSectors] = -1;
		      PzParityValues2[NbrMomentumSectors] = 0;
		      SzParityValues2[NbrMomentumSectors] = 0;
		    }
		  else
		    {
		      SzParityValues1[0] = 0;
		    }
		}
	      else
		{
		  // spin and valley
		  NbrSymmetrySectors = 1;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      NbrSymmetrySectors *= 2;
		    }
		  if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		    {
		      NbrSymmetrySectors *= 2;
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  // NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz + MinPz + MinEz);
		  // NbrParticlesBand1UpMinus[0] = (NbrParticles + MinSz - MinPz - MinEz);
		  // NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz + MinPz - MinEz);
		  // NbrParticlesBand1DownMinus[0] = (NbrParticles - MinSz - MinPz + MinEz);			  
		  NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz + MinPz);
		  NbrParticlesBand1UpMinus[0] = (NbrParticles + MinSz - MinPz);
		  NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz + MinPz);
		  NbrParticlesBand1DownMinus[0] = (NbrParticles - MinSz - MinPz);			  
		  // if ((NbrParticlesBand1UpPlus[0] < 0) || (NbrParticlesBand1UpMinus[0] < 0) || (NbrParticlesBand1DownPlus[0] < 0) || (NbrParticlesBand1DownMinus[0] < 0)
		  //     || ((NbrParticlesBand1UpPlus[0] & 3) != 0) ||  ((NbrParticlesBand1UpMinus[0] & 3) != 0)
		  //     || ((NbrParticlesBand1DownPlus[0] & 3) != 0) ||  ((NbrParticlesBand1DownMinus[0] & 3) != 0))
		  //   {
		  //     cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to 4N_{up,+}=" << NbrParticlesBand1UpPlus[0]
		  // 	   << " 4N_{up,-}=" << NbrParticlesBand1UpMinus[0] << " 4N_{down,+}=" << NbrParticlesBand1DownPlus[0]
		  // 	   << " 4N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		  //     return 0;
		  //   }
		  NbrParticlesBand1UpPlus[0] /= 4;
		  NbrParticlesBand1UpMinus[0] /= 4;
		  NbrParticlesBand1DownPlus[0] /= 4;
		  NbrParticlesBand1DownMinus[0] /= 4;
		  PzParityValues1[0] = 0;
		  SzParityValues1[0] = 0;
		  PzParityValues2[0] = 0;
		  SzParityValues2[0] = 0;
		  // if ((NbrParticlesBand1UpPlus[0] > NbrParticles) || (NbrParticlesBand1UpMinus[0] > NbrParticles)
		  //     || (NbrParticlesBand1DownPlus[0] > NbrParticles) || (NbrParticlesBand1DownMinus[0] > NbrParticles))
		  //   {
		  //     cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to N_{up,+}=" << NbrParticlesBand1UpPlus[0]
		  // 	   << " N_{up,-}=" << NbrParticlesBand1UpMinus[0] << " N_{down,+}=" << NbrParticlesBand1DownPlus[0]
		  // 	   << " N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		  //     return 0;
		  //   }
		  // cout << "N_{up,+}=" << NbrParticlesBand1UpPlus[0] << " N_{up,-}=" << NbrParticlesBand1UpMinus[0]
		  //      << " N_{down,+}=" << NbrParticlesBand1DownPlus[0] << " N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
			{
			  PzParityValues1[0] = 1;
			  SzParityValues1[0] = 1;
			  PzParityValues1[NbrMomentumSectors] = -1;
			  SzParityValues1[NbrMomentumSectors] = 1;
			  PzParityValues1[2 * NbrMomentumSectors] = 1;
			  SzParityValues1[2 * NbrMomentumSectors] = -1;
			  PzParityValues1[3 * NbrMomentumSectors] = -1;
			  SzParityValues1[3 * NbrMomentumSectors] = -1;
			  for (int i = NbrMomentumSectors; i < NbrSymmetrySectors; i += NbrMomentumSectors)
			    {
			      NbrParticlesBand1UpPlus[i] = NbrParticlesBand1UpPlus[0];
			      NbrParticlesBand1UpMinus[i] = NbrParticlesBand1UpMinus[0];
			      NbrParticlesBand1DownPlus[i] = NbrParticlesBand1DownPlus[0];
			      NbrParticlesBand1DownMinus[i] = NbrParticlesBand1DownMinus[0];
			      PzParityValues2[i] = 0;
			      SzParityValues2[i] = 0;
			    }
			}
		      else
			{
			  SzParityValues1[0] = 1;
			  NbrParticlesBand1UpPlus[NbrMomentumSectors] = NbrParticlesBand1UpPlus[0];
			  NbrParticlesBand1UpMinus[NbrMomentumSectors] = NbrParticlesBand1UpMinus[0];
			  NbrParticlesBand1DownPlus[NbrMomentumSectors] = NbrParticlesBand1DownPlus[0];
			  NbrParticlesBand1DownMinus[NbrMomentumSectors] = NbrParticlesBand1DownMinus[0];
			  PzParityValues1[NbrMomentumSectors] = 0;
			  SzParityValues1[NbrMomentumSectors] = -1;
			  PzParityValues2[NbrMomentumSectors] = 0;
			  SzParityValues2[NbrMomentumSectors] = 0;
			}
		    }
		  else
		    {
		      if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
			{
			  PzParityValues1[0] = 1;
			  NbrParticlesBand1UpPlus[NbrMomentumSectors] = NbrParticlesBand1UpPlus[0];
			  NbrParticlesBand1UpMinus[NbrMomentumSectors] = NbrParticlesBand1UpMinus[0];
			  NbrParticlesBand1DownPlus[NbrMomentumSectors] = NbrParticlesBand1DownPlus[0];
			  NbrParticlesBand1DownMinus[NbrMomentumSectors] = NbrParticlesBand1DownMinus[0];
			  PzParityValues1[NbrMomentumSectors] = -1;
			  SzParityValues1[NbrMomentumSectors] = 0;
			  PzParityValues2[NbrMomentumSectors] = 0;
			  SzParityValues2[NbrMomentumSectors] = 0;
			}
		    }
		}
	    }
	}
  
	KxMomenta = new int[NbrSymmetrySectors];
	KyMomenta = new int[NbrSymmetrySectors];
	SzValues = new int[NbrSymmetrySectors];
	PzValues = new int[NbrSymmetrySectors];
	EzValues = new int[NbrSymmetrySectors];
	int TmpIndex = 0;
	for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      for (int k = 0; k < NbrSymmetrySectors; k += NbrMomentumSectors)
		{
		  NbrParticlesBand1UpPlus[k + TmpIndex] = NbrParticlesBand1UpPlus[k];
		  NbrParticlesBand1UpMinus[k + TmpIndex] = NbrParticlesBand1UpMinus[k];
		  NbrParticlesBand1DownPlus[k + TmpIndex] = NbrParticlesBand1DownPlus[k];
		  NbrParticlesBand1DownMinus[k + TmpIndex] = NbrParticlesBand1DownMinus[k];
		  KxMomenta[k + TmpIndex] = i;
		  KyMomenta[k + TmpIndex] = j;
		  SzValues[k + TmpIndex] = MinSz;
		  PzValues[k + TmpIndex] = MinPz;
		  //		  EzValues[k + TmpIndex] = MinEz;
		  EzValues[k + TmpIndex] = 0;
		  SzParityValues1[k + TmpIndex] = SzParityValues1[k];
		  SzParityValues2[k + TmpIndex] = SzParityValues2[k];
		  PzParityValues1[k + TmpIndex] = PzParityValues1[k];
		  PzParityValues2[k + TmpIndex] = PzParityValues2[k];
		}
	      ++TmpIndex;
	    }
	}
    }
  else
    {
      // quantum number selection from an external file       
      MultiColumnASCIIFile SymmetrySectorsFile;
      int MininumNumberColumns = 0;
      if (SymmetrySectorsFile.Parse(Manager.GetString("selected-sectors")) == false)
	{
	  SymmetrySectorsFile.DumpErrors(cout) << endl;
	  return 0;
	}
      if (SymmetrySectorsFile.GetNbrLines() == 0)
	{
	  cout << Manager.GetString("selected-sectors") << " is an empty file" << endl;
	  return 0;
	}
      //      if (Manager.GetBoolean("conserve-bandoccuption") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  MininumNumberColumns = 2;
		  if (SymmetrySectorsFile.GetNbrColumns() < 2)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least two)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  if (Manager.GetBoolean("open-y") == true)
		    {
		      KyMomenta = new int [NbrSymmetrySectors];
		      for (int i = 0; i < NbrSymmetrySectors; ++i)
			{
			  KyMomenta[i] = 0;
			}
		    }
		  else
		    {
		      KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		    }
		  SzValues = new int [NbrSymmetrySectors];
		  PzValues = new int [NbrSymmetrySectors];
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      SzValues[i] = NbrParticles;
		      PzValues[i] = NbrParticles;
		      EzValues[i] = NbrParticles;
		      NbrParticlesBand1UpPlus[i] = NbrParticles;
		      NbrParticlesBand1UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = 0;
		      PzParityValues1[i] = 0;
		      SzParityValues1[i] = 0;
		      PzParityValues2[i] = 0;
		      SzParityValues2[i] = 0;
		    }
		}
	      else
		{
		  // valley, no spin
		  MininumNumberColumns = 3;
		  if (SymmetrySectorsFile.GetNbrColumns() < 3)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least three when using --add-valley)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  if (Manager.GetBoolean("open-y") == true)
		    {
		      KyMomenta = new int [NbrSymmetrySectors];
		      for (int i = 0; i < NbrSymmetrySectors; ++i)
			{
			  KyMomenta[i] = 0;
			}
		      PzValues = SymmetrySectorsFile.GetAsIntegerArray(1);
		    }
		  else
		    {
		      KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		      PzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		    }
		  SzValues = new int [NbrSymmetrySectors];
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      SzValues[i] = NbrParticles;
		      EzValues[i] = NbrParticles;
		      NbrParticlesBand1UpPlus[i] = (NbrParticles + PzValues[i]) / 2;
		      NbrParticlesBand1UpMinus[i] = (NbrParticles - PzValues[i]) / 2;
		      NbrParticlesBand1DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = 0;
		      PzParityValues1[i] = 0;
		      SzParityValues1[i] = 0;
		      PzParityValues2[i] = 0;
		      SzParityValues2[i] = 0;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  // spin, no valley
		  MininumNumberColumns = 3;
		  if (SymmetrySectorsFile.GetNbrColumns() < 3)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least three when using --add-spin)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  if (Manager.GetBoolean("open-y") == true)
		    {
		      KyMomenta = new int [NbrSymmetrySectors];
		      for (int i = 0; i < NbrSymmetrySectors; ++i)
			{
			  KyMomenta[i] = 0;
			}
		      SzValues = SymmetrySectorsFile.GetAsIntegerArray(1);
		    }
		  else
		    {
		      KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		      SzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		    }
		  PzValues = new int [NbrSymmetrySectors];
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      PzValues[i] = NbrParticles;
		      EzValues[i] = NbrParticles;
		      NbrParticlesBand1UpPlus[i] = (NbrParticles + SzValues[i]) / 2;
		      NbrParticlesBand1UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = (NbrParticles - PzValues[i]) / 2;
		      NbrParticlesBand1DownMinus[i] = 0;
		      PzParityValues1[i] = 0;
		      SzParityValues1[i] = 0;
		      PzParityValues2[i] = 0;
		      SzParityValues2[i] = 0;
		    }
		}
	      else
		{
		  // spin and valley
		  MininumNumberColumns = 5;
		  if (SymmetrySectorsFile.GetNbrColumns() < 5)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least five when using --add-spin and --add-valley)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  if (Manager.GetBoolean("open-y") == true)
		    {
		      KyMomenta = new int [NbrSymmetrySectors];
		      for (int i = 0; i < NbrSymmetrySectors; ++i)
			{
			  KyMomenta[i] = 0;
			}
		      SzValues = SymmetrySectorsFile.GetAsIntegerArray(1);
		      PzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		    }
		  else
		    {
		      KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		      SzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		      PzValues = SymmetrySectorsFile.GetAsIntegerArray(3);
		    }
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      NbrParticlesBand1UpPlus[i] = (NbrParticles + SzValues[i] + PzValues[i] + EzValues[i]) / 4;
		      NbrParticlesBand1UpMinus[i] = (NbrParticles + SzValues[i] - PzValues[i] - EzValues[i]) / 4;
		      NbrParticlesBand1DownPlus[i] = (NbrParticles - SzValues[i] + PzValues[i] - EzValues[i]) / 4;
		      NbrParticlesBand1DownMinus[i] = (NbrParticles - SzValues[i] - PzValues[i] + EzValues[i]) / 4;
		      EzValues[i] = 0;
		      PzParityValues1[i] = 0;
		      SzParityValues1[i] = 0;
		      PzParityValues2[i] = 0;
		      SzParityValues2[i] = 0;
		    }
		}
	    }
	}
      if (Manager.GetBoolean("variable-nbreigenvalues") == true)
	{
	  if (SymmetrySectorsFile.GetNbrColumns() < (MininumNumberColumns + 1))
	    {
	      cout << "--variable-nbreigenvalues requires at least " << (MininumNumberColumns + 1) << " columns in " << Manager.GetString("selected-sectors") << " (has only " << SymmetrySectorsFile.GetNbrColumns() << " columns)" << endl;
	      return 0;
	    }
	  NbrRequestedEigenstates = SymmetrySectorsFile.GetAsIntegerArray(MininumNumberColumns);
	}
    }
  
  int TotalDim = 0;
  for (int SymmetrySectorIndex = 0; SymmetrySectorIndex < NbrSymmetrySectors; ++SymmetrySectorIndex)
    {
      if (Manager.GetBoolean("add-valley") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("open-y") == false)
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ") : " << endl;
		}
	      else
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ") : " << endl;
		}
	    }
	  else
	    {
	      if (UseSzMinusSzSymmetry == false)
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",Sz<->-Sz=" << SzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",Sz<->-Sz=" << SzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (UsePzMinusPzSymmetry == false)
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",Pz<->-Pz=" << PzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",Pz<->-Pz=" << PzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	    }
	  else
	    {
	      if ((UsePzMinusPzSymmetry == false) && (UseSzMinusSzSymmetry == false))
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ") : " << endl;		      
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ") : " << endl;		      
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",Sz<->-Sz=" << SzParityValues1[SymmetrySectorIndex] << ",Pz<->-Pz=" << PzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",Sz<->-Sz=" << SzParityValues1[SymmetrySectorIndex] << ",Pz<->-Pz=" << PzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	    }
	}
      ParticleOnSphereWithSpin* Space = 0;
      AbstractQHEHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (Manager.GetBoolean("boson") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  ParticleOnSphere* TmpSpace = 0;
		  // no valley, no spin
		  if ((NbrSitesX * NbrSitesY) <= 64)
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  TmpSpace = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], Manager.GetString("hilbert-directory"));
			}
		      else
			{
			  TmpSpace = new FermionOnSquareOpenLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], Manager.GetString("hilbert-directory"));
			}
		      Space = new ParticleOnSphereWithPolarizedSpin(TmpSpace);
		    }
		  else
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  TmpSpace = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			}
		      else
			{
			  cout << "open bc in the y direction and long Hilbert space not available" << endl;
			  return 0;
			  //			  TmpSpace = new FermionOnOpenSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			}
		      Space = new ParticleOnSphereWithPolarizedSpin(TmpSpace);
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == true)
		    {
		      cout << "open bc in the y direction not available" << endl;
		      return 0;
		    }
		  // valley but no spin
		  if ((NbrSitesX * NbrSitesY) <= 32)
		    {
		      if (PzParityValues1[SymmetrySectorIndex] == 0)
			{
			  Space = new FermionOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (NbrParticles + PzValues[SymmetrySectorIndex]) / 2, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
										      10000000ul);
			}
		      else
			{
			  cout << "pz parity is not implemented" << endl; 
			  return 0;
			  // Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (NbrParticles, NbrSitesX, NbrSitesY,
			  // 									    KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
			  // 									    PzValues[SymmetrySectorIndex],
			  // 									    (PzParityValues1[SymmetrySectorIndex] == -1), 10000000ul);
			  // Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry (NbrParticles, NbrSitesX, NbrSitesY,
			  // 											  KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
			  // 											  PzValues[SymmetrySectorIndex],
			  // 											  (PzParityValues1[SymmetrySectorIndex] == -1), 10000000ul);
			}			    
		    }
		  else
		    {
		      Space = new FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong (NbrParticles, (NbrParticles + PzValues[SymmetrySectorIndex]) / 2, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], 10000000ul);
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("open-y") == true)
		{
		  cout << "open bc in the y direction not available" << endl;
		  return 0;
		}
	      // spinful case
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if ((NbrSitesX * NbrSitesY) <= 32)
		    {
		      if (SzParityValues1[SymmetrySectorIndex] == 0)
			{
			  Space = new FermionOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (NbrParticles + SzValues[SymmetrySectorIndex]) / 2, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
										      10000000ul);
			}
		      else
			{
			  cout << "sz parity is not implemented" << endl;
			  return 0;
			  // Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (NbrParticles, NbrSitesX, NbrSitesY,
			  // 									    KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
			  // 									    SzValues[SymmetrySectorIndex],
			  // 									    (SzParityValues1[SymmetrySectorIndex] == -1), 10000000ul);
			}
		    }
		  else
		    {
		      Space = new FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong (NbrParticles, (NbrParticles + SzValues[SymmetrySectorIndex]) / 2, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], 10000000ul);
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("open-y") == true)
		    {
		      cout << "open bc in the y direction not available" << endl;
		      return 0;
		    }
		  if ((NbrSitesX * NbrSitesY) <= 16)
		    {
		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
										  SzValues[SymmetrySectorIndex], PzValues[SymmetrySectorIndex], 10000000ul);
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 32)
			{			  
			  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
											  SzValues[SymmetrySectorIndex], PzValues[SymmetrySectorIndex], 10000000ul);
			}
		      else
			{
			  cout << "SU(4) not supported with more than 32 momenta" << endl;
			  Space = 0;			      
			}
		    }
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("open-y") == true)
	    {
	      cout << "open bc in the y direction not available" << endl;
	      return 0;
	    }
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  // no valley, no spin
		  ParticleOnSphere* TmpSpace;
		  if (((NbrSitesX * NbrSitesY) + (NbrParticles - 1)) <= 64)
		    {
		      TmpSpace = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		    }
		  else
		    {
		      TmpSpace = new BosonOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		    }
		  Space = new ParticleOnSphereWithPolarizedSpin(TmpSpace);
		}
	      else
		{
		  // valley but no spin
		  cout << "bosons are not supported" << endl;
		  // Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		}
	    }
	  else
	    {
	      cout << "bosons are not supported" << endl;
	      return 0;
	    }
	}
      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
      TotalDim += Space->GetHilbertSpaceDimension();
      if ((Space->GetHilbertSpaceDimension() > 0)  && (Manager.GetBoolean("hilbert-only") == false))
	{
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  // if (Manager.GetBoolean("real-interaction"))
	  //   {
	  //     if (Manager.GetBoolean("add-valley") == false)
	  // 	{
	  // 	  Hamiltonian = new ParticleOnLatticeFromFileInteractionThreeBandRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
	  // 											Manager.GetString("interaction-file"),
	  // 											TightBindingModel, Manager.GetBoolean("flat-band"), 
	  // 											Manager.GetDouble("interaction-rescaling"),
	  // 											Manager.GetBoolean("add-spin"),
	  // 											Architecture.GetArchitecture(), Memory);
	  // 	}
	  //     else
	  // 	{
	  // 	  Hamiltonian = new ParticleOnLatticeFromFileInteractionThreeBandWithSpinRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
	  // 												Manager.GetString("interaction-file"),
	  // 												TightBindingModel, Manager.GetBoolean("flat-band"), 
	  // 												Manager.GetDouble("interaction-rescaling"),
	  // 												Manager.GetBoolean("add-spin"),
	  // 												Architecture.GetArchitecture(), Memory);
	  // 	}		
	  //   }
	  // else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if (Manager.GetString("threebody-interaction-file") ==0)
		    {
		      Hamiltonian = new ParticleOnLatticeFromFileInteractionOneBandHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
												Manager.GetString("interaction-file"),
												TightBindingModel, Manager.GetBoolean("flat-band"), 
												Manager.GetDouble("interaction-rescaling"),
												Manager.GetBoolean("add-spin"),
												Manager.GetBoolean("open-y"),
												Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      Hamiltonian = new ParticleOnLatticeFromFileInteractionOneBandThreeBodyHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
													 Manager.GetString("threebody-interaction-file"), Manager.GetDouble("threebody-interaction-rescaling"),
													 Manager.GetString("interaction-file"), Manager.GetDouble("interaction-rescaling"),
													 TightBindingModel, Manager.GetBoolean("flat-band"), 
													 
													 Manager.GetBoolean("add-spin"),
													 Architecture.GetArchitecture(), Memory);
		    }
		}
	      else
		{
		  if (Manager.GetString("threebody-interaction-file") ==0)
		    {
		      Hamiltonian = new ParticleOnLatticeFromFileInteractionOneBandWithSpinHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
													Manager.GetString("interaction-file"),
													TightBindingModel, Manager.GetBoolean("flat-band"), 
													Manager.GetDouble("interaction-rescaling"),
													Manager.GetBoolean("add-spin"),
													Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      cout << "three-body interaction for spinful case is not yet implemented" << endl;
		      return 0;
		    }
		}				
	    }
	  
	  
	  char* ContentPrefix = new char[256];
	  char* EigenstateOutputFile = new char [512];
	  char* TmpExtention = new char[256];
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  if (Manager.GetBoolean("open-y") == false)
		    {
		      sprintf (ContentPrefix, "%d %d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d_ky_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d", KxMomenta[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d", KxMomenta[SymmetrySectorIndex]);
		    }
		}
	      else
		{
		  if (UseSzMinusSzSymmetry == false)
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  sprintf (ContentPrefix, "%d %d %d", SzValues[SymmetrySectorIndex], KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d", SzValues[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_sz_%d", KxMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  sprintf (ContentPrefix, "%d %d %d %d", SzValues[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d_szsym_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   SzValues[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d", SzValues[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_sz_%d_szsym_%d", KxMomenta[SymmetrySectorIndex],
				   SzValues[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex]);
			}
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  if (UsePzMinusPzSymmetry == false)
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  sprintf (ContentPrefix, "%d %d %d", PzValues[SymmetrySectorIndex], KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], PzValues[SymmetrySectorIndex]);
			}
		      else
			{
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  sprintf (ContentPrefix, "%d %d %d %d", PzValues[SymmetrySectorIndex], PzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_pzsym_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], PzParityValues1[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d", PzValues[SymmetrySectorIndex], PzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_pz_%d_pzsym_%d", KxMomenta[SymmetrySectorIndex], 
				   PzValues[SymmetrySectorIndex], PzParityValues1[SymmetrySectorIndex]);
			}
		    }
		}
	      else
		{
		  if ((UsePzMinusPzSymmetry == false) && (UseSzMinusSzSymmetry == false))
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  sprintf (ContentPrefix, "%d %d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_sz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_pz_%d_sz_%d", KxMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("open-y") == false)
			{
			  sprintf (ContentPrefix, "%d %d %d %d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   PzParityValues1[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_sz_%d_pzsym_%d_szsym_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   PzParityValues1[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   PzParityValues1[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_pz_%d_sz_%d_pzsym_%d_szsym_%d", KxMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   PzParityValues1[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex]);
			}
		    }
		}
	    }
	  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	  delete[] TmpExtention;
	  Hamiltonian->ShiftHamiltonian(EnergyShift);
	  // if (Manager.GetBoolean("real-interaction"))
	  //   {
	  //     GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix,
	  // 			       CommentLine, EnergyShift,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  //    if (NbrRequestedEigenstates != 0)
	  //	{
	  //	  Task.SetNbrEigenvalues(NbrRequestedEigenstates[SymmetrySectorIndex]);
	  //	}
	  //     MainTaskOperation TaskOperation (&Task);
	  //     TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  //   }
	  // else
	    {
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix,
					  CommentLine, EnergyShift,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      if (NbrRequestedEigenstates != 0)
		{
		  Task.SetNbrEigenvalues(NbrRequestedEigenstates[SymmetrySectorIndex]);
		}
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	  FirstRunFlag = false;
	  
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}
      delete Space;
    }
  cout << "Total dim=" << TotalDim << endl;
  return 0;
}

