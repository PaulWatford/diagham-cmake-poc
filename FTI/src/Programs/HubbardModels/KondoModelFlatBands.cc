#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeWithSpinAndMagneticImpurityRealSpace.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"


#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
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
  cout.precision(14);
  OptionManager Manager ("KondoModelFlatBands" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 2);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-orbitals", "number of orbitals", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-value", "twice the spin Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption ('\n', "testhermitian-error", "error threshold when testing hermiticy (0 for machine accuracy)", 0.0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type KondoModelFlatBands -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrOrbitals = Manager.GetInteger("nbr-orbitals"); 
  int TotalSz = Manager.GetInteger("sz-value");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (SzSymmetryFlag == false)
    {
      sprintf (StatisticPrefix, "fermions_kondoflatband");
    }
  else
    {
      sprintf (StatisticPrefix, "fermions_kondoflatband_szsym");
    }

  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_x_%d_%d_n_%d", StatisticPrefix, NbrOrbitals, NbrParticles);
  
  char* CommentLine = new char [256];
  if (SzSymmetryFlag == false)
    {
      sprintf (CommentLine, "sz");
    }
  else
    {
      sprintf (CommentLine, "sz szp");
    }

  char* EigenvalueOutputFile = new char [512];
  sprintf(EigenvalueOutputFile, "%s_sz_%d.dat", FilePrefix, TotalSz);


  int SzParitySector = -1;
  int MaxSzParitySector = 1;
  if (SzSymmetryFlag == false)
    {
      SzParitySector = 1;
    }
  else
    {
      if ((Manager.GetInteger("sz-parity") != 0) && (TotalSz == 0))
	{
	  SzParitySector = Manager.GetInteger("sz-parity");
	  MaxSzParitySector = SzParitySector;
	}
    }

  bool FirstRunFlag = true;

  for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
    {
      ParticleOnSphereWithSpin* Space = 0;
      AbstractHamiltonian* Hamiltonian = 0;
      if (SzSymmetryFlag == false)
	{
	  cout << "Sz = " << TotalSz <<  endl;
	  Space = new FermionOnLatticeWithSpinAndMagneticImpurityRealSpace (NbrParticles, TotalSz, NbrOrbitals);
	}
      else
	{
	  bool MinusParitySector = true;
	  if (SzParitySector == 1)
	    MinusParitySector = false;
	  cout << "Sz = " << TotalSz << "  SzParity = " << SzParitySector<< endl;
	  //	  Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, TotalSz, NbrSites, MinusParitySector);
	}
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      
      // Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceHamiltonian(Space, NbrParticles + 1, NbrOrbitals + 1,
      // 								      TightBindingMatrix, TightBindingMatrix,
      // 								      DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
      // 								      DensityDensityInteractionupdown, 
      // 								      Architecture.GetArchitecture(), Memory);
      
      char* ContentPrefix = new char[256];
      if (SzSymmetryFlag == false)
	{
	  sprintf (ContentPrefix, "%d", TotalSz);
	}
      else
	{
	  sprintf (ContentPrefix, "%d %d", TotalSz, SzParitySector);
	}
      char* EigenstateOutputFile;
      char* TmpExtention = new char [512];
      if (SzSymmetryFlag == false)
	{
	  sprintf (TmpExtention, "_sz_%d", TotalSz);
	}
      else
	{
	  sprintf (TmpExtention, "_szp_%d_sz_%d", SzParitySector, TotalSz);
	}
      char* TmpExtentionSpectrum = new char [32];
      sprintf (TmpExtentionSpectrum, "_sz_%d.dat", TotalSz);
      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, TmpExtentionSpectrum, TmpExtention);
      
      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
      FirstRunFlag = false;
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      cout << "------------------------------------" << endl;
      delete Hamiltonian;
      delete Space;
      delete[] EigenstateOutputFile;
      delete[] ContentPrefix;
    }

  return 0;
}
