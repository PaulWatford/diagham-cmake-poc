#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzGutzwillerProjection.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"

#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 
#include <limits>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(std::numeric_limits<double>::max_digits10);

  // some running options and help
  OptionManager Manager ("FQHESphereTimeEvolutionTwoBodyGeneric" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, false, 2);
  ParticleManager.AddOptionGroup(&Manager);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption('\n', "initial-state", "name of the file containing the initial vector upon which e^{-iHt} acts");
  (*SystemGroup) += new BooleanOption  ('\n', "complex", "initial vector is a complex vector");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-energy", "compute the energy of each time-evolved vector");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('A', "all-sz", "assume a hilbert space including all sz values");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "no-spin", "do not compute the S^2 value of the state");
  (*SystemGroup) += new BooleanOption  ('\n', "no-szparity", "do not compute the parity under the Sz<->-Sz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
//  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
//  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "show-extracted", "show values extracted from file name");

  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "sma");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potentials");
 
  (*SystemGroup) += new SingleDoubleOption ('\n', "time-step", "time interval between two snap shots", 0.1);
  (*SystemGroup) += new SingleDoubleOption ('\n', "time-shift", "time shift to add in the vectors name, if initial state for current run is not at t = 0", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-steps", "number of points to evaluate", 1);
  
  (*SystemGroup) += new SingleDoubleOption ('\n', "precision", "convergence precision", 1.0e-14);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximal number of iterations", 100);
  
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereTimeEvolutionTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  int PairParity = Manager.GetInteger("pair-parity");
  bool SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
  bool SzMinusParity = Manager.GetBoolean("minus-szparity");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool LzMinusParity = Manager.GetBoolean("minus-lzparity");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  int TmpTotalSz = TotalSz;
  bool AllSzFlag = false;
  bool GutzwillerFlag = false;
  if (Manager.GetBoolean("all-sz"))
    {
      TmpTotalSz = -1;
      AllSzFlag = true;
    }
  int LzSymmetry = 0;
  int SzSymmetry = 0;

  if (NbrParticles == 0)
    {
      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, TmpTotalSz, SzSymmetry, LzSymmetry, FermionFlag, AllSzFlag) == false)
	{
	  return -1;
	}
      else
	{
	  if (strstr(Manager.GetString("state"), "_gutzwiller_"))
	    {
	      GutzwillerFlag = true;
	    }
	  if (SzSymmetry != 0)
	    {
	      SzSymmetrizedBasis = true;
	      if (SzSymmetry < 0)
		{
		  SzMinusParity = true;
		}
	      else
		{
		  SzMinusParity = false;
		}
	    }
	  if (LzSymmetry != 0)
	    {
	      LzSymmetrizedBasis = true;
	      if (LzSymmetry < 0)
		{
		  LzMinusParity = true;
		}
	      else
		{
		  LzMinusParity = false;
		}
	    }
	    //	  if (!Manager.GetBoolean("all-sz"))
	    //	    TotalSz=TmpTotalSz;
	  if (Manager.GetBoolean("show-extracted") == true)
	    {
	      cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz;
	      if (AllSzFlag == false)
		{		  
		  cout<< "  TotalSz=" << TotalSz;
		}
	      else
		{
		  cout<< "  All Sz sectors ";
		}
	      if (LzSymmetrizedBasis == true)
		{
		  cout << "  Lz symmetrized basis ";
		  if (LzMinusParity == true)
		    cout << "(minus parity) ";
		  else
		    cout << "(plus parity) ";
		}
	      if (SzSymmetrizedBasis == true)
		{
		  cout << "  Sz symmetrized basis ";
		  if (SzMinusParity == true)
		    cout << "(minus parity) ";
		else
		  cout << "(plus parity) ";
		}
	      cout << endl;
	    }
	}
    }
  if (Manager.GetBoolean("lzsymmetrized-basis") == true)
    {
      LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
      LzMinusParity = Manager.GetBoolean("minus-lzparity");      
    }
  if (Manager.GetBoolean("szsymmetrized-basis") == true)
    {
      SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
      SzMinusParity = Manager.GetBoolean("minus-szparity");
    }
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	{
	  if ((strcmp ("bosons", Manager.GetString("statistics")) == 0))
	    {
	      FermionFlag = false;
	    }
	  else
	    {
	      cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	    }  
	}
    }
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  
  double TmpTime = Manager.GetDouble("time-step");
  int NbrTimeSteps = Manager.GetInteger("nbr-steps");
  double TimeShift = Manager.GetDouble("time-shift");
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("initial-state"),
						   NbrParticles, LzMax, TotalLz, FermionFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("initial-state") << endl;
      return -1;
    }
  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << LzMax << " Lz=" << TotalLz << " ";

  bool onDiskCacheFlag = Manager.GetBoolean("allow-disk-storage");
  double** PseudoPotentials  = new double*[10];
  for (int i = 0; i < 3; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
	PseudoPotentials[i][j] = 0.0;
    };
  double** OneBodyPseudoPotentials  = new double*[3];
  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double * OneBodyPotentialUpDown = 0;

  int NbrUp = (NbrParticles + TotalSz) >> 1;
  int NbrDown = (NbrParticles - TotalSz) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  if (Manager.GetBoolean("all-sz") == true)
   {
     NbrUp = NbrParticles;
     NbrDown = 0;
   } 

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereSU2GetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials,
					   OneBodyPseudoPotentials) == false)
	return -1;
    }

  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
 
  char* OutputNamePrefix = new char [512];
  char* NormName = new char[512];
  char* EnergyName;
  char* InteractionName = new char[strlen(Manager.GetString("interaction-name")) + 1];
  strcpy(InteractionName, Manager.GetString("interaction-name"));
  
  ParticleOnSphereWithSpin* Space = (ParticleOnSphereWithSpin*)ParticleManager.GetHilbertSpace(TotalLz);

  if (FermionFlag == false)
    {
      sprintf (OutputNamePrefix, "bosons_sphere_su2_%s_n_%d_2s_%d_t", Manager.GetString("interaction-name"), NbrParticles, LzMax);
      sprintf (NormName, "bosons_sphere_su2_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_norm.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      if (Manager.GetBoolean("compute-energy"))
      {
	EnergyName = new char[512];
	sprintf (EnergyName, "bosons_sphere_su2_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_energy.dat", InteractionName, NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      }
      
    }
  else
    {
      sprintf (OutputNamePrefix, "fermions_sphere_su2_%s_n_%d_2s_%d_t", Manager.GetString("interaction-name"), NbrParticles, LzMax);
      sprintf (NormName, "fermions_sphere_su2_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_norm.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      if (Manager.GetBoolean("compute-energy"))
      {
	EnergyName = new char[512];
	sprintf (EnergyName, "fermions_sphere_su2_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_energy.dat", InteractionName, NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      }
    }

  delete[] InteractionName;
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  char* StateFileName = Manager.GetString("initial-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  ComplexVector TmpInitialState (Space->GetHilbertSpaceDimension());
  if (Manager.GetBoolean("complex") == false)
  {
    RealVector InputState;
    if (InputState.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
    if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< Space->GetHilbertSpaceDimension() << endl;
      return -1;
    }
    TmpInitialState = InputState;
  }
  else
    {
    ComplexVector InputState;
    if (InputState.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
    if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< Space->GetHilbertSpaceDimension() << endl;
      return -1;
    }
    TmpInitialState = InputState;
  }
  
  
    
  cout << " Initial state Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
  
  
  
  ParticleOnSphereWithSpinGenericHamiltonian* Hamiltonian = new ParticleOnSphereWithSpinGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
								   OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], OneBodyPseudoPotentials[2], 
								   Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);
											 
//    double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
	
  ofstream File;
  File.open(NormName, ios::binary | ios::out);
  File.precision(std::numeric_limits<double>::max_digits10);
  File << "# t Norm dNorm" << endl;
  
  ofstream FileEnergy;
  if (Manager.GetBoolean("compute-energy"))
  {
    FileEnergy.open(EnergyName, ios::binary | ios::out);
    FileEnergy.precision(std::numeric_limits<double>::max_digits10);
    FileEnergy << "# t E "<< endl;
  }
  
  double Norm;
  int TmpExpansionOrder;
  ComplexVector TmpState (Space->GetHilbertSpaceDimension()) ;
  ComplexVector TmpState1 (Space->GetHilbertSpaceDimension()) ;
  Complex TmpCoefficient;
  for (int i = 0; i < NbrTimeSteps; ++i)
  {
    TmpState.Copy(TmpInitialState);
    Norm = TmpState.Norm();
    double TmpNorm = 1.0;
    TmpExpansionOrder = 0;
    TmpCoefficient = 1.0;
    cout << "Computing state " << (i + 1) << "/" << NbrTimeSteps << " at t = " << (TmpTime * i) << endl;
    while (((fabs(TmpNorm) > Manager.GetDouble("precision")) || (TmpExpansionOrder < 1)) && (TmpExpansionOrder <= Manager.GetInteger("iter-max")))
    {
      TmpExpansionOrder += 1;
      TmpCoefficient = -TmpCoefficient * TmpTime * Complex(0.0, 1.0) / ((double) TmpExpansionOrder);
      VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpState), (&TmpState1));
      Operation.ApplyOperation(Architecture.GetArchitecture());
      TmpState.Copy(TmpState1);
      TmpNorm = sqrt(TmpCoefficient.Re*TmpCoefficient.Re + TmpCoefficient.Im*TmpCoefficient.Im) * TmpState.Norm();
      AddComplexLinearCombinationOperation Operation1 (&TmpInitialState, &TmpState1, 1, &TmpCoefficient);
      Operation1.ApplyOperation(Architecture.GetArchitecture());
      Norm = TmpInitialState.Norm();
      
      cout << "Norm = " << Norm << " +/- " << TmpNorm << " for step " << TmpExpansionOrder << endl;      
    }
    File << (TimeShift + (i + 1)*TmpTime) << " " << Norm << " " << TmpNorm << endl;
    cout << endl;  
    
    char* OutputName = new char [strlen(OutputNamePrefix)+ 16];
    sprintf (OutputName, "%s_%g_lz_%d.vec", OutputNamePrefix, TimeShift + (i + 1)*TmpTime, TotalLz);
    TmpInitialState.WriteVector(OutputName);
    delete[] OutputName;
    
    if (Manager.GetBoolean("compute-energy"))
    {
      VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpInitialState), (&TmpState));
      Operation.ApplyOperation(Architecture.GetArchitecture());
      Complex Energy = (TmpInitialState) * (TmpState);
      FileEnergy << (TimeShift + (i + 1)*TmpTime) << " " << Energy.Re << endl;
      cout << "E = " << Energy.Re << endl;
    }
  }
  
  File.close();
  if (Manager.GetBoolean("compute-energy"))
  {
    FileEnergy.close();
    delete[] EnergyName;
  }
  delete Hamiltonian;
  delete[] OutputNamePrefix;
  delete[] NormName;
  
 
  return 0;
}
