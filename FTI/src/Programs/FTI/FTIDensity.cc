#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpaceLong.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpaceLong.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU6SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU12SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong.h"

#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSquareLatticeWithGenericSpinBandDensityOperator.h"
#include "Operator/ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator.h"

#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// extract the system information from the file name
//
// stateFileName = state file name
// nbrBands = number of bands
// flag3d, flagDecoupled, flagWannier = values of flags
// nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz, totalSpin, statistics = references on system parameters
// maxBand0, maxBand1, maxBand2, maxBand3, minBand0, minBand1, minBand2, minBand3 = references on system parameters
// return value = true if no error occured
bool FTIDensityGetSystemInformationFromState(char* stateFileName, bool flag3d, bool flagDecoupled, bool flagWannier, int nbrBands,
					     int& nbrParticles, int& nbrSitesX, int& nbrSitesY, int& nbrSiteZ, int& totalKx, int& totalKy, int& totalKz, int& totalSpin, bool& Statistics,
					     int& maxBand0, int& maxBand1, int& maxBand2, int& maxBand3, int& minBand0, int& minBand1, int& minBand2, int& minBand3);

// get the Hilbert for any of the potential cases
//
// flag3d, flagDecoupled, flagWannier = values of flags
// nbrBands = number of bands
// nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz, totalSpin, statistics = system parameters
// maxBand0, maxBand1, maxBand2, maxBand3, minBand0, minBand1, minBand2, minBand3 = system parameters
// manager = reference on the option manager
// return value = pointer to the Hilbert space (null if an error occured)
 
ParticleOnSphere* FTIDensityGetHilbertSpace(bool flag3d, bool flagDecoupled, bool flagWannier, int nbrBands, int nbrParticles, int nbrSitesX, int nbrSitesY, int nbrSiteZ,
					    int totalKx, int totalKy, int totalKz, int totalSpin, bool statistics,
					    int maxBand0, int maxBand1, int maxBand2, int maxBand3, int minBand0, int minBand1, int minBand2, int minBand3,
					    OptionManager& manager);

// compute the momentum creation index with momentum transfer from the right state to the left state
//
// creationMomentumIndex = linearized
// leftTotalKx = momentum along the x direction for the left state
// leftTotalKy = momentum along the y direction for the left state
// leftTotalKz = momentum along the z direction for the left state
// totalKx = momentum along the x direction for the right state
// totalKy = momentum along the y direction for the right state
// totalKz = momentum along the z direction for the right state
// creationKx = reference on the momentum along the x direction for creation operator
// creationKy = reference on the momentum along the y direction for creation operator
// creationKz = reference on the momentum along the z direction for creation operator
// NbrSitesX = number of unit cells along the x direction
// NbrSitesY = number of unit cells along the y direction
// NbrSitesZ = number of unit cells along the z direction
// return value =momentum creation index including momentum transfer
int FTIDensityComputeMomentumTransfer(int creationMomentumIndex, int leftTotalKx, int leftTotalKy, int leftTotalKz,
				      int totalKx, int totalKy, int totalKz, int& creationKx, int& creationKy, int& creationKz,
				      int nbrSitesX, int nbrSitesY, int nbrSitesZ);


int main(int argc, char** argv)
{
  OptionManager Manager ("FTIDensity" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleStringOption  ('\n', "left-states", "single column file describing a degenerated ground state for the bra (aka laft states of the expectation values), should be used with --off-diagonal");
  (*SystemGroup) += new BooleanOption  ('\n', "off-diagonal", "use --degenerated-groundstate as a list of states and compute all the cross density terms");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume that the FTI states are made of two decoupled FCI copies");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sigma", "in decoupled mode, only compute in a given spin/sigma sector if non-negative", -1);
  (*SystemGroup) += new BooleanOption ('\n', "time-reversal", "apply complex conjugation to the left state before computing expectation values");
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption  ('\n', "Wannier", "Wannier basis");
  (*SystemGroup) += new BooleanOption  ('\n', "rhorho", "also compute the density-density expectation values");
  (*SystemGroup) += new BooleanOption  ('\n', "intraband-only", "when computing the density-density expectation values, only consider the intra band terms");
  (*SystemGroup) += new SingleStringOption ('\n', "allowed-orbitals", "provide an ASCII file indicating which orbitals are allowed");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rho.dat extension");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIDensity -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  
  int NbrRightStates = 1;
  int NbrLeftStates = 1;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  ComplexVector* RightGroundStates = 0;
  char** RightGroundStateFiles = 0;
  ComplexVector* LeftGroundStates = 0;
  char** LeftGroundStateFiles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int* TotalKz = 0;
  int* LeftTotalKx = 0;
  int* LeftTotalKy = 0;
  int* LeftTotalKz = 0;
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  int NbrSitesZ = 0;
  int MaxBand0 = -1;
  int MaxBand1 = -1;
  int MaxBand2 = -1;
  int MaxBand3 = -1;
  int MinBand0 = 0;
  int MinBand1 = 0;
  int MinBand2 = 0;
  int MinBand3 = 0;
  int SelectSigmaSector = Manager.GetInteger("sigma");
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  bool Flag3d = Manager.GetBoolean("3d");
  bool FlagDecoupled = Manager.GetBoolean("decoupled");
  int TotalSpin = 0;
  int NbrBands = Manager.GetInteger("nbr-subbands");
  bool FlagWannier = Manager.GetBoolean("Wannier");
  

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FTIDensity -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      RightGroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      TotalKz = new int[1];
      Coefficients = new double[1];
      RightGroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (RightGroundStateFiles[0], Manager.GetString("ground-file"));
      Coefficients[0] = 1.0;
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrRightStates = DegeneratedFile.GetNbrLines();
       RightGroundStateFiles = new char* [NbrRightStates];
       TotalKx = new int[NbrRightStates];
       TotalKy = new int[NbrRightStates];
       TotalKz = new int[NbrRightStates];
       for (int i = 0; i < NbrRightStates; ++i)
	 {
	   RightGroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (RightGroundStateFiles[i], DegeneratedFile(0, i));		   
	 }
       if (DegeneratedFile.GetNbrColumns() == 1)
	 {
	   Coefficients = new double[NbrRightStates];
	   for (int i = 0; i < NbrRightStates; ++i)
	     Coefficients[i] = 1.0 / ((double) NbrRightStates);
	 }
       else
	 {
	   double TmpSum = 0.0;
	   Coefficients = DegeneratedFile.GetAsDoubleArray(1);
	   for (int i = 0; i < NbrRightStates; ++i)
	     TmpSum += Coefficients[i];
	   TmpSum = 1.0 / TmpSum;
	   for (int i = 0; i < NbrRightStates; ++i)
	     Coefficients[i] *= TmpSum;
	 }
    }


  if (Manager.GetString("left-states") != 0)
   {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("left-states")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrLeftStates = DegeneratedFile.GetNbrLines();
       LeftGroundStateFiles = new char* [NbrLeftStates];
       LeftTotalKx = new int[NbrLeftStates];
       LeftTotalKy = new int[NbrLeftStates];
       LeftTotalKz = new int[NbrLeftStates];
       for (int i = 0; i < NbrLeftStates; ++i)
	 {
	   LeftGroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (LeftGroundStateFiles[i], DegeneratedFile(0, i));		   
	 }
    }



  
  for (int i = 0; i < NbrRightStates; ++i)
    {
      if (FTIDensityGetSystemInformationFromState(RightGroundStateFiles[i], Flag3d, FlagDecoupled, FlagWannier, NbrBands,
						  NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, TotalKx[i], TotalKy[i], TotalKz[i], TotalSpin, Statistics,
						  MaxBand0, MaxBand1, MaxBand2, MaxBand3, MinBand0, MinBand1, MinBand2, MinBand3) == false)
	{
	  cout << "error while retrieving system parameters from file name " << RightGroundStateFiles[i] << endl;
	  return -1;
	}
    }

  if (Manager.GetString("left-states") != 0)
    {
      for (int i = 0; i < NbrRightStates; ++i)
	{
	  int TmpNbrParticles = 0;
	  int TmpNbrSitesX = 0, TmpNbrSitesY = 0, TmpNbrSitesZ = 0, TmpTotalSpin = 0;
	  int TmpMaxBand0 = 0, TmpMaxBand1 = 0, TmpMaxBand2 = 0, TmpMaxBand3 = 0;
	  int TmpMinBand0 = 0, TmpMinBand1 = 0, TmpMinBand2 = 0, TmpMinBand3 = 0;
	  bool TmpStatistics = false;
	  if (FTIDensityGetSystemInformationFromState(LeftGroundStateFiles[i], Flag3d, FlagDecoupled, FlagWannier, NbrBands,
						      TmpNbrParticles, TmpNbrSitesX, TmpNbrSitesY, TmpNbrSitesZ, LeftTotalKx[i], LeftTotalKy[i], LeftTotalKz[i], TmpTotalSpin, TmpStatistics,
						      TmpMaxBand0, TmpMaxBand1, TmpMaxBand2, TmpMaxBand3, TmpMinBand0, TmpMinBand1, TmpMinBand2, TmpMinBand3) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << RightGroundStateFiles[i] << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpNbrSitesX != NbrSitesX) || (TmpNbrSitesY != NbrSitesY) || (TmpNbrSitesZ != NbrSitesZ) || (TmpTotalSpin != TotalSpin)
	      || (TmpMaxBand0 != MaxBand0) || (TmpMaxBand1 != MaxBand1) || (TmpMaxBand2 != MaxBand2) || (TmpMaxBand3 != MaxBand3)
	      || (TmpMinBand0 != MinBand0) || (TmpMinBand1 != MinBand1) || (TmpMinBand2 != MinBand2) || (TmpMinBand3 != MinBand3)
	      || (TmpStatistics != Statistics))
	    {
	      cout << "left and right states have different system parameters: " << endl;
	      cout << "Statistics=";
	      if (TmpStatistics == true)
		{
		  cout << "fermions";
		}
	      else
		{
		  cout << "bosons";
		}
	      cout << " vs ";
	      if (Statistics == true)
		{
		  cout << "fermions";
		}
	      else
		{
		  cout << "bosons";
		}
	      cout << endl;
	      cout << "NbrParticles=" << TmpNbrParticles << " vs " << NbrParticles << ", NbrSitesX=" << TmpNbrSitesX << " vs " << NbrSitesX
		   << ", NbrSitesY=" << TmpNbrSitesY << " vs " << NbrSitesY << ", NbrSitesZ=" << TmpNbrSitesZ << " vs " << NbrSitesZ
		   << ", TotalSpin=" << TmpTotalSpin << " vs " << TotalSpin
		   << ", MaxBand0=" << TmpMaxBand0 << " vs " << MaxBand0 << ", MaxBand1=" << TmpMaxBand1 << " vs " << MaxBand1 << ", MaxBand2=" << TmpMaxBand2 << " vs " << MaxBand2 << ", MaxBand3=" << TmpMaxBand3 << " vs " << MaxBand3
		   << ", MinBand0=" << TmpMinBand0 << " vs " << MinBand0 << ", MinBand1=" << TmpMinBand1 << " vs " << MinBand1 << ", MinBand2=" << TmpMinBand2 << " vs " << MinBand2 << ", MinBand3=" << TmpMinBand3 << " vs " << MinBand3 << endl;
	      return -1;
	    }
	}
    }
  
  RightGroundStates = new ComplexVector [NbrRightStates];  
  int TotalNbrSites;
  if(FlagWannier == false || (FlagWannier == true &&  TotalKx[0]>-1))
    TotalNbrSites = NbrSitesX * NbrSitesY * NbrSitesZ;
  else
    TotalNbrSites = NbrSitesY * NbrSitesZ;
  int* NbrRightGroundStatePerMomentumSector = new int[TotalNbrSites];
  int* NbrLeftGroundStatePerMomentumSector = new int[TotalNbrSites];
  int* NbrGroundStatePerMomentumSector = new int[TotalNbrSites];
  ComplexVector** RightGroundStatePerMomentumSector = new ComplexVector*[TotalNbrSites];
  ComplexVector** LeftGroundStatePerMomentumSector = new ComplexVector*[TotalNbrSites];
  double** CoefficientPerMomentumSector = new double*[TotalNbrSites];
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      NbrGroundStatePerMomentumSector[i] = 0;
      NbrRightGroundStatePerMomentumSector[i] = 0;
      RightGroundStatePerMomentumSector[i] = 0;
      NbrLeftGroundStatePerMomentumSector[i] = 0;
      LeftGroundStatePerMomentumSector[i] = 0;
      CoefficientPerMomentumSector[i] = 0;
    }
  for (int i = 0; i < NbrRightStates; ++i)
    {
      if (RightGroundStates[i].ReadVector (RightGroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << RightGroundStateFiles[i] << endl;
	  return -1;      
	}
      int TmpIndex;
      if ((FlagWannier == false) || ((FlagWannier == true) && (TotalKx[0]>-1)))
	TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSitesZ + TotalKz[i];
      NbrRightGroundStatePerMomentumSector[TmpIndex]++; 
    }
  if (Manager.GetString("left-states") != 0)
    {
      for (int i = 0; i < NbrLeftStates; ++i)
	{
	  if (LeftGroundStates[i].ReadVector (LeftGroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << LeftGroundStateFiles[i] << endl;
	      return -1;      
	    }
	  int TmpIndex;
	  if ((FlagWannier == false) || ((FlagWannier == true) && (TotalKx[0]>-1)))
	    TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
	  else
	    TmpIndex = TotalKy[i] * NbrSitesZ + TotalKz[i];
	  NbrLeftGroundStatePerMomentumSector[TmpIndex]++; 
	}
    }

  
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      if (NbrRightGroundStatePerMomentumSector[i] > 0)
	{
	  RightGroundStatePerMomentumSector[i] = new ComplexVector[NbrRightGroundStatePerMomentumSector[i]];
	  CoefficientPerMomentumSector[i] = new double[NbrRightGroundStatePerMomentumSector[i]];
	}
      NbrRightGroundStatePerMomentumSector[i] = 0;
      if (NbrLeftGroundStatePerMomentumSector[i] > 0)
	{
	  LeftGroundStatePerMomentumSector[i] = new ComplexVector[NbrLeftGroundStatePerMomentumSector[i]];
	  CoefficientPerMomentumSector[i] = new double[NbrLeftGroundStatePerMomentumSector[i]];
	}
      NbrLeftGroundStatePerMomentumSector[i] = 0;
    }
  
  for (int i = 0; i < NbrRightStates; ++i)
    {
      int TmpIndex;
      if ((FlagWannier == false) || ((FlagWannier == true) && (TotalKx[0]>-1)))
	TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSitesZ + TotalKz[i];
      RightGroundStatePerMomentumSector[TmpIndex][NbrRightGroundStatePerMomentumSector[TmpIndex]] = RightGroundStates[i];
      CoefficientPerMomentumSector[TmpIndex][NbrRightGroundStatePerMomentumSector[TmpIndex]] = Coefficients[i];
      NbrRightGroundStatePerMomentumSector[TmpIndex]++;
    }  
  if (Manager.GetString("left-states") != 0)
    {
      for (int i = 0; i < NbrLeftStates; ++i)
	{
	  int TmpIndex;
	  if ((FlagWannier == false) || ((FlagWannier == true) && (TotalKx[0]>-1)))
	    TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
	  else
	    TmpIndex = TotalKy[i] * NbrSitesZ + TotalKz[i];
	  LeftGroundStatePerMomentumSector[TmpIndex][NbrLeftGroundStatePerMomentumSector[TmpIndex]] = LeftGroundStates[i];
	  NbrLeftGroundStatePerMomentumSector[TmpIndex]++;
	}
    }

  int NbrDensityIndices = 0;
  int NbrDensityPartialTraces = 0;
  int* CreationMomentumIndices = 0;
  int* AnnihilationMomentumIndices = 0;
  int* CreationSigmaIndices = 0;
  int* AnnihilationSigmaIndices = 0;
  char** IndexLabels = 0;
  char** PartialTraceLabels = 0;
  char* FileHeader = new char[256];
  
  int MaxNbrSpaces; 
  if ((FlagWannier == false)  || ((FlagWannier == true) && (TotalKx[0] > -1)))
    MaxNbrSpaces = NbrSitesX * NbrSitesY * NbrSitesZ;
  else
    MaxNbrSpaces = NbrSitesY * NbrSitesZ;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrRightStates; ++i)
    {
      int TmpIndex;
      if ((FlagWannier == false)  || ((FlagWannier == true) && (TotalKx[0] > -1)))
	TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSitesZ + TotalKz[i];
      if (Spaces[TmpIndex] == 0)
	{
	  Spaces[TmpIndex] = FTIDensityGetHilbertSpace(Flag3d, FlagDecoupled, FlagWannier, NbrBands,
						       NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, TotalKx[i], TotalKy[i], TotalKz[i], TotalSpin, Statistics,
						       MaxBand0, MaxBand1, MaxBand2, MaxBand3, MinBand0, MinBand1, MinBand2, MinBand3, Manager);
	  if (Spaces[TmpIndex] == 0)
	    {
	      return 0;
	    }
	  if (Flag3d == false)
	    {
	      if (NbrBands == 1)
		{
		  sprintf (FileHeader, "# kx ky <c^+ c>");
		  NbrDensityIndices = NbrSitesX * NbrSitesY;
		  NbrDensityPartialTraces = 1;
		  PartialTraceLabels = new char* [NbrDensityPartialTraces];
		  for (int i = 0; i < NbrDensityPartialTraces; ++i)
		    {
		      PartialTraceLabels[i] = new char [256];
		      sprintf(PartialTraceLabels[i], "");
		    }
		  CreationMomentumIndices = new int[NbrDensityIndices];
		  AnnihilationMomentumIndices = new int[NbrDensityIndices];
		  CreationSigmaIndices = new int[NbrDensityIndices];
		  AnnihilationSigmaIndices = new int[NbrDensityIndices];
		  IndexLabels = new char*[NbrDensityIndices];
		  NbrDensityIndices = 0;
		  for (int kx = 0; kx < NbrSitesX; ++kx)
		    {	
		      for (int ky = 0; ky < NbrSitesY; ++ky)
			{
			  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
			  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
			  CreationSigmaIndices[NbrDensityIndices] = 0;
			  AnnihilationSigmaIndices[NbrDensityIndices] = 0;
			  IndexLabels[NbrDensityIndices] = new char[256];
			  sprintf(IndexLabels[NbrDensityIndices], "%d %d", kx, ky);
			  ++NbrDensityIndices;
			}
		    }
		}
	      if (NbrBands == 2)
		{
		  if (FlagDecoupled == false)
		    {
		      sprintf (FileHeader, "# kx ky sigma sigma' <c^+ c>");
		      NbrDensityIndices = 4 * NbrSitesX * NbrSitesY;
		      NbrDensityPartialTraces = 2;
		      PartialTraceLabels = new char* [NbrDensityPartialTraces];
		      for (int i = 0; i < NbrDensityPartialTraces; ++i)
			{
			  PartialTraceLabels[i] = new char [256];
			  sprintf(PartialTraceLabels[i], "sigma=%d", i);
			}
		      CreationMomentumIndices = new int[NbrDensityIndices];
		      AnnihilationMomentumIndices = new int[NbrDensityIndices];
		      CreationSigmaIndices = new int[NbrDensityIndices];
		      AnnihilationSigmaIndices = new int[NbrDensityIndices];
		      IndexLabels = new char*[NbrDensityIndices];
		      NbrDensityIndices = 0;
		      for (int kx = 0; kx < NbrSitesX; ++kx)
			{	
			  for (int ky = 0; ky < NbrSitesY; ++ky)
			    {
			      for (int i = 0; i <= 1; ++i)
				{
				  for (int j = 0; j <= 1; ++j)
				    {
				      CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				      AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				      CreationSigmaIndices[NbrDensityIndices] = i;
				      AnnihilationSigmaIndices[NbrDensityIndices] = j;
				      IndexLabels[NbrDensityIndices] = new char[256];
				      sprintf(IndexLabels[NbrDensityIndices], "%d %d %d %d", kx, ky, i, j);
				      ++NbrDensityIndices;
				    }
				}
			    }
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("off-diagonal") == false)
			{
			  sprintf (FileHeader, "# kx ky sigma <c^+ c>");			  
			}
		      else
			{
			  sprintf (FileHeader, "# psi_i phi_j kx1 ky1 kx2 ky2 sigma <psi_i | c^+_{kx1,ky1} c_{kx2,ky2} | phi_j>");
			}
		      int MinSigma = 0;
		      int MaxSigma = 1;		      
		      if ((SelectSigmaSector >= 0) && (SelectSigmaSector < 2))
			{
			  MinSigma = SelectSigmaSector;
			  MaxSigma = SelectSigmaSector;
			}
		      NbrDensityIndices = 2 * NbrSitesX * NbrSitesY;
		      NbrDensityPartialTraces = 2;
		      PartialTraceLabels = new char* [NbrDensityPartialTraces];
		      for (int i = 0; i < NbrDensityPartialTraces; ++i)
			{
			  PartialTraceLabels[i] = new char [256];
			  sprintf(PartialTraceLabels[i], "sigma=%d", MinSigma + i);
			}
		      CreationMomentumIndices = new int[NbrDensityIndices];
		      AnnihilationMomentumIndices = new int[NbrDensityIndices];
		      CreationSigmaIndices = new int[NbrDensityIndices];
		      AnnihilationSigmaIndices = new int[NbrDensityIndices];
		      IndexLabels = new char*[NbrDensityIndices];
		      NbrDensityIndices = 0;
		      for (int kx = 0; kx < NbrSitesX; ++kx)
			{	
			  for (int ky = 0; ky < NbrSitesY; ++ky)
			    {
			      for (int i = MinSigma; i <= MaxSigma; ++i)
				{
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = i;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d %d", kx, ky, i);
				  ++NbrDensityIndices;
				}
			    }
			}
		    }
		}
	      if (NbrBands == 3)
		{
		  if (FlagDecoupled == false)
		    {
		      sprintf (FileHeader, "# kx ky sigma sigma' <c^+ c>");
		      NbrDensityIndices = 9 * NbrSitesX * NbrSitesY;
		      NbrDensityPartialTraces = 3;
		      PartialTraceLabels = new char* [NbrDensityPartialTraces];
		      for (int i = 0; i < NbrDensityPartialTraces; ++i)
			{
			  PartialTraceLabels[i] = new char [256];
			  sprintf(PartialTraceLabels[i], "sigma=%d", i);
			}
		      CreationMomentumIndices = new int[NbrDensityIndices];
		      AnnihilationMomentumIndices = new int[NbrDensityIndices];
		      CreationSigmaIndices = new int[NbrDensityIndices];
		      AnnihilationSigmaIndices = new int[NbrDensityIndices];
		      IndexLabels = new char*[NbrDensityIndices];
		      NbrDensityIndices = 0;
		      for (int kx = 0; kx < NbrSitesX; ++kx)
			{	
			  for (int ky = 0; ky < NbrSitesY; ++ky)
			    {
			      for (int i = 0; i <= 2; ++i)
				{
				  for (int j = 0; j <= 2; ++j)
				    {
				      CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				      AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				      CreationSigmaIndices[NbrDensityIndices] = i;
				      AnnihilationSigmaIndices[NbrDensityIndices] = j;
				      IndexLabels[NbrDensityIndices] = new char[256];
				      sprintf(IndexLabels[NbrDensityIndices], "%d %d %d %d", kx, ky, i, j);
				      ++NbrDensityIndices;
				    }
				}
			    }
			}
		    }
		}
	      if (NbrBands == 4)
		{
		  sprintf (FileHeader, "# kx ky spin sigma sigma' <c^+ c>");
		  NbrDensityIndices = 8 * NbrSitesX * NbrSitesY;
		  NbrDensityPartialTraces = 4;
		  PartialTraceLabels = new char* [NbrDensityPartialTraces];
		  for (int i = 0; i < NbrDensityPartialTraces; i++)
		    {
		      PartialTraceLabels[i] = new char [256];
		      sprintf(PartialTraceLabels[i], "spin=%d sigma=%d", (i >> 1), (i & 1));
		    }
		  CreationMomentumIndices = new int[NbrDensityIndices];
		  AnnihilationMomentumIndices = new int[NbrDensityIndices];
		  CreationSigmaIndices = new int[NbrDensityIndices];
		  AnnihilationSigmaIndices = new int[NbrDensityIndices];
		  IndexLabels = new char*[NbrDensityIndices];
		  NbrDensityIndices = 0;
		  for (int kx = 0; kx < NbrSitesX; ++kx)
		    {	
		      for (int ky = 0; ky < NbrSitesY; ++ky)
			{
			  for (int i = 0; i <= 1; ++i)
			    {
			      for (int j = 0; j <= 1; ++j)
				{
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = j;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d 0 %d %d", kx, ky, i, j);
				  ++NbrDensityIndices;
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = 2 + i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = 2 + j;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d 1 %d %d", kx, ky, i, j);
				  ++NbrDensityIndices;
				}
			    }
			}
		    }
		}
	      if (NbrBands == 6)
		{
		  sprintf (FileHeader, "# kx ky spin sigma sigma' <c^+ c>");
		  NbrDensityIndices = 18 * NbrSitesX * NbrSitesY;
		  NbrDensityPartialTraces = 6;
		  PartialTraceLabels = new char* [NbrDensityPartialTraces];
		  for (int i = 0; i < NbrDensityPartialTraces; i++)
		    {
		      PartialTraceLabels[i] = new char [256];
		      sprintf(PartialTraceLabels[i], "spin=%d sigma=%d", (i / 3), (i % 3));
		    }
		  CreationMomentumIndices = new int[NbrDensityIndices];
		  AnnihilationMomentumIndices = new int[NbrDensityIndices];
		  CreationSigmaIndices = new int[NbrDensityIndices];
		  AnnihilationSigmaIndices = new int[NbrDensityIndices];
		  IndexLabels = new char*[NbrDensityIndices];
		  NbrDensityIndices = 0;
		  for (int kx = 0; kx < NbrSitesX; ++kx)
		    {	
		      for (int ky = 0; ky < NbrSitesY; ++ky)
			{
			  for (int i = 0; i <= 2; ++i)
			    {
			      for (int j = 0; j <= 2; ++j)
				{
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = j;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d 0 %d %d", kx, ky, i, j);
				  ++NbrDensityIndices;
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSitesY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = 3 + i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = 3 + j;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d 1 %d %d", kx, ky, i, j);
				  ++NbrDensityIndices;
				}
			    }
			}
		    }
		}	      
	    }
	  else
	    {
	    }
	}
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != RightGroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space (" << Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() << ") and ground state (" << RightGroundStates[i].GetLargeVectorDimension() << ")" << endl;
	  return 0;
	}
    }
  
  if (Manager.GetString("left-states") != 0)
    {
      for (int i = 0; i < NbrLeftStates; ++i)
	{
	  int TmpIndex;
	  if ((FlagWannier == false)  || ((FlagWannier == true) && (LeftTotalKx[0] > -1)))
	    TmpIndex = (((LeftTotalKx[i] * NbrSitesY) + LeftTotalKy[i]) * NbrSitesZ) + LeftTotalKz[i];
	  else
	    TmpIndex = LeftTotalKy[i] * NbrSitesZ + LeftTotalKz[i];
	  if (Spaces[TmpIndex] == 0)
	    {
	      Spaces[TmpIndex] = FTIDensityGetHilbertSpace(Flag3d, FlagDecoupled, FlagWannier, NbrBands,
							   NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, LeftTotalKx[i], LeftTotalKy[i], LeftTotalKz[i], TotalSpin, Statistics,
							   MaxBand0, MaxBand1, MaxBand2, MaxBand3, MinBand0, MinBand1, MinBand2, MinBand3, Manager);
	      if (Spaces[TmpIndex] == 0)
		{
		  return 0;
		}
	    }
	}
    }
  
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
    }
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(RightGroundStateFiles[0], "vec", "rho.dat");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << RightGroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  File << FileHeader << endl;

  if (Manager.GetBoolean("off-diagonal") == false)
    {
      if (Flag3d == false)
	{
	  if (NbrBands == 1)
	    {
	      for (int i = 0; i < NbrRightStates; ++i)
		{
		  Complex TmpTotalDensity = 0.0;
		  Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		  for (int j = 0; j < NbrDensityPartialTraces; ++j)
		    {
		      PartialTraces[j] = 0.0;
		    }
		  int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
		  for (int j = 0 ; j < NbrDensityIndices; ++j)
		    {
		      ParticleOnSphereDensityOperator TmpOperator ((ParticleOnSphere*) Spaces[TmpIndex], CreationMomentumIndices[j], AnnihilationMomentumIndices[j]);
		      Complex TmpElement = TmpOperator.MatrixElement(RightGroundStates[i], RightGroundStates[i]);
		      if (CreationSigmaIndices[j] == AnnihilationSigmaIndices[j])
			{
			  PartialTraces[CreationSigmaIndices[j]] += TmpElement;
			  TmpTotalDensity += TmpElement;
			}
		      File << IndexLabels[j] << " " << TmpElement << endl;
		    }
		  for (int j = 0; j < NbrDensityPartialTraces; ++j)
		    {
		      File << "# partial density " << PartialTraceLabels[j] << " = " << PartialTraces[j] << endl;
		    }
		  File << "# total density = " << TmpTotalDensity << endl;
		}
	    }
	  else
	    {
	      if (NbrBands >= 2)
		{
		  for (int i = 0; i < NbrRightStates; ++i)
		    {
		      Complex TmpTotalDensity = 0.0;
		      Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		      for (int j = 0; j < NbrDensityPartialTraces; ++j)
			{
			  PartialTraces[j] = 0.0;
			}
		      int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
		      for (int j = 0 ; j < NbrDensityIndices; ++j)
			{
			  ParticleOnSquareLatticeWithGenericSpinBandDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], CreationMomentumIndices[j], CreationSigmaIndices[j], AnnihilationMomentumIndices[j], AnnihilationSigmaIndices[j]);
			  Complex TmpElement = TmpOperator.MatrixElement(RightGroundStates[i], RightGroundStates[i]);
			  if (CreationSigmaIndices[j] == AnnihilationSigmaIndices[j])
			    {
			      PartialTraces[CreationSigmaIndices[j]] += TmpElement;
			      TmpTotalDensity += TmpElement;
			    }
			  File << IndexLabels[j] << " " << TmpElement << endl;
			}
		      for (int j = 0; j < NbrDensityPartialTraces; ++j)
			{
			  File << "# partial density " << PartialTraceLabels[j] << " = " << PartialTraces[j] << endl;
			}
		      File << "# total density = " << TmpTotalDensity << endl;
		    }
		}
	    }
	}
    }
  else
    {
      // off diagonal case
      if (Flag3d == false)
	{
	  if (NbrBands == 1)
	    {
	      for (int i = 0; i < NbrRightStates; ++i)
		{
		}
	      cout << "Error, --off-diagonal is not implemented for a single band, no spin/valley Hilbert space" << endl;
	      return 0;
	    }
	  else
	    {
	      if (NbrBands >= 2)
		{
		  if (Manager.GetString("left-states") == 0)
		    {
		      for (int i = 0; i < NbrRightStates; ++i)
			{
			  Complex TmpTotalDensity = 0.0;
			  Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
			  for (int j = 0; j < NbrDensityPartialTraces; ++j)
			    {
			      PartialTraces[j] = 0.0;
			    }
			  int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
			  for (int j = 0; j < NbrRightStates; ++j)
			    {
			      int TmpLeftIndex = (((TotalKx[j] * NbrSitesY) + TotalKy[j]) * NbrSitesZ) + TotalKz[j];
			      Spaces[TmpIndex]->SetTargetSpace(Spaces[TmpLeftIndex]);
			      for (int k = 0 ; k < NbrDensityIndices; ++k)
				{
				  int TmpKx;
				  int TmpKy;
				  int TmpKz;
				  int TmpCreationMomentumIndex = FTIDensityComputeMomentumTransfer(CreationMomentumIndices[k], TotalKx[j], TotalKy[j], TotalKz[j],
												   TotalKx[i], TotalKy[i], TotalKz[i],
												   TmpKx, TmpKy, TmpKz,
												   NbrSitesX, NbrSitesY, NbrSitesZ);
				  ParticleOnSquareLatticeWithGenericSpinBandDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], TmpCreationMomentumIndex,
													 CreationSigmaIndices[k], AnnihilationMomentumIndices[k],
													 AnnihilationSigmaIndices[k]);
				  Complex TmpElement;
				  if (Manager.GetBoolean("time-reversal") == false)
				    {
				      OperatorMatrixElementOperation TmpOperation (&TmpOperator, RightGroundStates[j], RightGroundStates[i]);
				      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
				      TmpElement = TmpOperation.GetScalar();
					//				      TmpElement = TmpOperator.MatrixElement(RightGroundStates[j], RightGroundStates[i]);
				    }
				  else
				    {
				      TmpElement = TmpOperator.ConjugateMatrixElement(RightGroundStates[j], RightGroundStates[i]);
				    }
				  if ((CreationSigmaIndices[k] == AnnihilationSigmaIndices[k]) && (TmpCreationMomentumIndex == AnnihilationMomentumIndices[k]) && (i == j))
				    {
				      PartialTraces[CreationSigmaIndices[k]] += TmpElement;
				      TmpTotalDensity += TmpElement;
				    }
				  if (Flag3d == false)
				    {
				      File << j << " " << i << " " << TmpKx << " " << TmpKy << " " << IndexLabels[k] << " " << TmpElement << endl;
				    }
				  else
				    {
				      File << j << " " << i << " " << TmpKx << " " << TmpKy << " " << TmpKz << " " << IndexLabels[k] << " " << TmpElement << endl;
				    }
				}
			    }
			  for (int j = 0; j < NbrDensityPartialTraces; ++j)
			    {
			      File << "# partial density " << PartialTraceLabels[j] << " = " << PartialTraces[j] << endl;
			    }
			  File << "# total density = " << TmpTotalDensity << endl;
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrRightStates; ++i)
			{
			  int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
			  for (int j = 0; j < NbrLeftStates; ++j)
			    {
			      int TmpLeftIndex = (((LeftTotalKx[j] * NbrSitesY) + LeftTotalKy[j]) * NbrSitesZ) + LeftTotalKz[j];
			      Spaces[TmpIndex]->SetTargetSpace(Spaces[TmpLeftIndex]);
			      for (int k = 0 ; k < NbrDensityIndices; ++k)
				{
				  int TmpKx;
				  int TmpKy;
				  int TmpKz;
				  int TmpCreationMomentumIndex = FTIDensityComputeMomentumTransfer(CreationMomentumIndices[k], LeftTotalKx[j], LeftTotalKy[j], LeftTotalKz[j],
												   TotalKx[i], TotalKy[i], TotalKz[i],
												   TmpKx, TmpKy, TmpKz,
												   NbrSitesX, NbrSitesY, NbrSitesZ);
				  ParticleOnSquareLatticeWithGenericSpinBandDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], TmpCreationMomentumIndex, CreationSigmaIndices[k], AnnihilationMomentumIndices[k], AnnihilationSigmaIndices[k]);
				  Complex TmpElement;
				  if (Manager.GetBoolean("time-reversal") == false)
				    {
				      TmpElement = TmpOperator.MatrixElement(LeftGroundStates[j], RightGroundStates[i]);
				    }
				  else
				    {
				      TmpElement = TmpOperator.ConjugateMatrixElement(LeftGroundStates[j], RightGroundStates[i]);
				    }
				  if (TmpIndex != TmpLeftIndex)
				    {				  
				      if (Flag3d == false)
					{
					  File << j << " " << i << TmpKx << " " << TmpKy << " " << IndexLabels[k] << " " << TmpElement << endl;
					}
				      else
					{
					  File << j << " " << i << TmpKx << " " << TmpKy << " " << TmpKz << " " << IndexLabels[k] << " " << TmpElement << endl;
					}
				    }
				  else
				    {
				      File << j << " " << i << " " << IndexLabels[k] << " " << TmpElement << endl;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  File.close();

  if (Manager.GetBoolean("rhorho") == true)
    {
      int NbrDensityDensityIndices = 0;
      int* CreationMomentumIndices1 = 0;
      int* CreationSigmaIndices1 = 0;
      int* CreationMomentumIndices2 = 0;
      int* CreationSigmaIndices2 = 0;
      int* AnnihilationMomentumIndices1 = 0;
      int* AnnihilationSigmaIndices1 = 0;
      int* AnnihilationMomentumIndices2 = 0;
      int* AnnihilationSigmaIndices2 = 0;
      char** DensityDensityIndexLabels = 0;
      char* DensityDensityFileHeader = new char[512];

      if (Flag3d == false)
	{
	  if (NbrBands == 1)
	    {
	    }
	  if (NbrBands == 2)
	    {
	      if (FlagDecoupled == false)
		{
		  sprintf (FileHeader, "# kx1 ky1 sigma1 kx2 ky2 sigma2 kx3 ky3 sigma3 kx4 ky4 sigma4 <c^+ c^+ c c>");
		  NbrDensityDensityIndices = 16 * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY);
		  CreationMomentumIndices1 = new int[NbrDensityDensityIndices];
		  CreationSigmaIndices1 = new int[NbrDensityDensityIndices];
		  CreationMomentumIndices2 = new int[NbrDensityDensityIndices];
		  CreationSigmaIndices2 = new int[NbrDensityDensityIndices];
		  AnnihilationMomentumIndices1 = new int[NbrDensityDensityIndices];
		  AnnihilationSigmaIndices1 = new int[NbrDensityDensityIndices];
		  AnnihilationMomentumIndices2 = new int[NbrDensityDensityIndices];
		  AnnihilationSigmaIndices2 = new int[NbrDensityDensityIndices];
		  IndexLabels = new char*[NbrDensityDensityIndices];
		  NbrDensityDensityIndices = 0;
		  for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
		    {	
		      for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
			{
			  for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
			    {	
			      for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
				{
				  for (int kx3 = 0; kx3 < NbrSitesX; ++kx3)
				    {	
				      for (int ky3 = 0; ky3 < NbrSitesY; ++ky3)
					{
					  int kx4 = kx1 + kx2 - kx3;
					  if (kx4 < 0)
					    {
					      kx4 += NbrSitesX;
					    }
					  kx4 %= NbrSitesX;
					  int ky4 = ky1 + ky2 - ky3;
					  if (ky4 < 0)
					    {
					      ky4 += NbrSitesY;
					    }
					  ky4 %= NbrSitesY;
					  for (int i = 0; i <= 1; ++i)
					    {
					      for (int j = 0; j <= 1; ++j)
						{
						  for (int k = 0; k <= 1; ++k)
						    {
						      for (int l = 0; l <= 1; ++l)
							{
							  CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
							  CreationSigmaIndices1[NbrDensityDensityIndices] = i;
							  CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
							  CreationSigmaIndices2[NbrDensityDensityIndices] = j;
							  AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
							  AnnihilationSigmaIndices1[NbrDensityDensityIndices] = k;
							  AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
							  AnnihilationSigmaIndices2[NbrDensityDensityIndices] = l;
							  IndexLabels[NbrDensityDensityIndices] = new char[512];
							  sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d %d %d %d %d %d %d %d %d %d %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
							  ++NbrDensityDensityIndices;
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	      else
		{
		  sprintf (FileHeader, "# kx1 ky1 sigma1 kx2 ky2 sigma2 kx3 ky3 sigma3 kx4 ky4 sigma4 <c^+ c^+ c c>");
		  NbrDensityDensityIndices = 3 * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY);
		  CreationMomentumIndices1 = new int[NbrDensityDensityIndices];
		  CreationSigmaIndices1 = new int[NbrDensityDensityIndices];
		  CreationMomentumIndices2 = new int[NbrDensityDensityIndices];
		  CreationSigmaIndices2 = new int[NbrDensityDensityIndices];
		  AnnihilationMomentumIndices1 = new int[NbrDensityDensityIndices];
		  AnnihilationSigmaIndices1 = new int[NbrDensityDensityIndices];
		  AnnihilationMomentumIndices2 = new int[NbrDensityDensityIndices];
		  AnnihilationSigmaIndices2 = new int[NbrDensityDensityIndices];
		  IndexLabels = new char*[NbrDensityDensityIndices];
		  NbrDensityDensityIndices = 0;
		  for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
		    {	
		      for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
			{
			  for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
			    {	
			      for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
				{
				  for (int kx3 = 0; kx3 < NbrSitesX; ++kx3)
				    {	
				      for (int ky3 = 0; ky3 < NbrSitesY; ++ky3)
					{
					  int kx4 = kx1 + kx2 - kx3;
					  if (kx4 < 0)
					    {
					      kx4 += NbrSitesX;
					    }
					  kx4 %= NbrSitesX;
					  int ky4 = ky1 + ky2 - ky3;
					  if (ky4 < 0)
					    {
					      ky4 += NbrSitesY;
					    }
					  ky4 %= NbrSitesY;
					  
					  CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
					  CreationSigmaIndices1[NbrDensityDensityIndices] = 0;
					  CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
					  CreationSigmaIndices2[NbrDensityDensityIndices] = 0;
					  AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
					  AnnihilationSigmaIndices1[NbrDensityDensityIndices] = 0;
					  AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
					  AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 0;
					  IndexLabels[NbrDensityDensityIndices] = new char[512];
					  sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d %d %d %d %d %d %d %d %d %d %d", kx1, ky1, 0, kx2, ky2, 0, kx3, ky3, 0, kx4, ky4, 0);
					  ++NbrDensityDensityIndices;

					  CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
					  CreationSigmaIndices1[NbrDensityDensityIndices] = 1;
					  CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
					  CreationSigmaIndices2[NbrDensityDensityIndices] = 1;
					  AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
					  AnnihilationSigmaIndices1[NbrDensityDensityIndices] = 1;
					  AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
					  AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 1;
					  IndexLabels[NbrDensityDensityIndices] = new char[512];
					  sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d %d %d %d %d %d %d %d %d %d %d", kx1, ky1, 1, kx2, ky2, 1, kx3, ky3, 1, kx4, ky4, 1);
					  ++NbrDensityDensityIndices;

					  CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
					  CreationSigmaIndices1[NbrDensityDensityIndices] = 0;
					  CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
					  CreationSigmaIndices2[NbrDensityDensityIndices] = 1;
					  AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
					  AnnihilationSigmaIndices1[NbrDensityDensityIndices] = 0;
					  AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
					  AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 1;
					  IndexLabels[NbrDensityDensityIndices] = new char[512];
					  sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d %d %d %d %d %d %d %d %d %d %d", kx1, ky1, 0, kx2, ky2, 1, kx3, ky3, 0, kx4, ky4, 1);
					  ++NbrDensityDensityIndices;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	  if (NbrBands == 3)
	    {
	      if (FlagDecoupled == false)
		{
		  if (Manager.GetBoolean("intraband-only") == true)
		    {		      
		      sprintf (FileHeader, "# kx1 ky1 sigma kx2 ky2 sigma kx3 ky3 sigma kx4 ky4 sigma <c^+ c^+ c c>");
		      NbrDensityDensityIndices = 3 * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY);
		      CreationMomentumIndices1 = new int[NbrDensityDensityIndices];
		      CreationSigmaIndices1 = new int[NbrDensityDensityIndices];
		      CreationMomentumIndices2 = new int[NbrDensityDensityIndices];
		      CreationSigmaIndices2 = new int[NbrDensityDensityIndices];
		      AnnihilationMomentumIndices1 = new int[NbrDensityDensityIndices];
		      AnnihilationSigmaIndices1 = new int[NbrDensityDensityIndices];
		      AnnihilationMomentumIndices2 = new int[NbrDensityDensityIndices];
		      AnnihilationSigmaIndices2 = new int[NbrDensityDensityIndices];
		      IndexLabels = new char*[NbrDensityDensityIndices];
		      NbrDensityDensityIndices = 0;
		      for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
			{	
			  for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
			    {
			      for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
				{	
				  for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
				    {
				      for (int kx3 = 0; kx3 < NbrSitesX; ++kx3)
					{	
					  for (int ky3 = 0; ky3 < NbrSitesY; ++ky3)
					    {
					      int kx4 = kx1 + kx2 - kx3;
					      if (kx4 < 0)
						{
						  kx4 += NbrSitesX;
						}
					      kx4 %= NbrSitesX;
					      int ky4 = ky1 + ky2 - ky3;
					      if (ky4 < 0)
						{
						  ky4 += NbrSitesY;
						}
					      ky4 %= NbrSitesY;
					      for (int i = 0; i <= 2; ++i)
						{
						  CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						  CreationSigmaIndices1[NbrDensityDensityIndices] = i;
						  CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						  CreationSigmaIndices2[NbrDensityDensityIndices] = i;
						  AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						  AnnihilationSigmaIndices1[NbrDensityDensityIndices] = i;
						  AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						  AnnihilationSigmaIndices2[NbrDensityDensityIndices] = i;
						  IndexLabels[NbrDensityDensityIndices] = new char[512];
						  sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d %d %d %d %d %d %d %d %d %d %d", kx1, ky1, i, kx2, ky2, i, kx3, ky3, i, kx4, ky4, i);
						  ++NbrDensityDensityIndices;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		  else
		    {
		      sprintf (FileHeader, "# kx1 ky1 sigma1 kx2 ky2 sigma2 kx3 ky3 sigma3 kx4 ky4 sigma4 <c^+ c^+ c c>");
		      NbrDensityDensityIndices = 81 * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY);
		      CreationMomentumIndices1 = new int[NbrDensityDensityIndices];
		      CreationSigmaIndices1 = new int[NbrDensityDensityIndices];
		      CreationMomentumIndices2 = new int[NbrDensityDensityIndices];
		      CreationSigmaIndices2 = new int[NbrDensityDensityIndices];
		      AnnihilationMomentumIndices1 = new int[NbrDensityDensityIndices];
		      AnnihilationSigmaIndices1 = new int[NbrDensityDensityIndices];
		      AnnihilationMomentumIndices2 = new int[NbrDensityDensityIndices];
		      AnnihilationSigmaIndices2 = new int[NbrDensityDensityIndices];
		      IndexLabels = new char*[NbrDensityDensityIndices];
		      NbrDensityDensityIndices = 0;
		      for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
			{	
			  for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
			    {
			      for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
				{	
				  for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
				    {
				      for (int kx3 = 0; kx3 < NbrSitesX; ++kx3)
					{	
					  for (int ky3 = 0; ky3 < NbrSitesY; ++ky3)
					    {
					      int kx4 = kx1 + kx2 - kx3;
					      if (kx4 < 0)
						{
						  kx4 += NbrSitesX;
						}
					      kx4 %= NbrSitesX;
					      int ky4 = ky1 + ky2 - ky3;
					      if (ky4 < 0)
						{
						  ky4 += NbrSitesY;
						}
					      ky4 %= NbrSitesY;
					      for (int i = 0; i <= 2; ++i)
						{
						  for (int j = 0; j <= 2; ++j)
						    {
						      for (int k = 0; k <= 2; ++k)
							{
							  for (int l = 0; l <= 2; ++l)
							    {
							      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
							      CreationSigmaIndices1[NbrDensityDensityIndices] = i;
							      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
							      CreationSigmaIndices2[NbrDensityDensityIndices] = j;
							      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
							      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = k;
							      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
							  AnnihilationSigmaIndices2[NbrDensityDensityIndices] = l;
							  IndexLabels[NbrDensityDensityIndices] = new char[512];
							  sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d %d %d %d %d %d %d %d %d %d %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
							  ++NbrDensityDensityIndices;
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	  if (NbrBands == 4)
	    {
	      sprintf (FileHeader, "# kx1 ky1 spin1 sigma1 kx2 ky2 spin2 sigma2 kx3 ky3 spin3 sigma3 kx4 ky4 spin4 sigma4 <c^+ c^+ c c>");
	      NbrDensityDensityIndices = 48* (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY);
	      CreationMomentumIndices1 = new int[NbrDensityDensityIndices];
	      CreationSigmaIndices1 = new int[NbrDensityDensityIndices];
	      CreationMomentumIndices2 = new int[NbrDensityDensityIndices];
	      CreationSigmaIndices2 = new int[NbrDensityDensityIndices];
	      AnnihilationMomentumIndices1 = new int[NbrDensityDensityIndices];
	      AnnihilationSigmaIndices1 = new int[NbrDensityDensityIndices];
	      AnnihilationMomentumIndices2 = new int[NbrDensityDensityIndices];
	      AnnihilationSigmaIndices2 = new int[NbrDensityDensityIndices];
	      IndexLabels = new char*[NbrDensityDensityIndices];
	      NbrDensityDensityIndices = 0;
	      
	      for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
		{	
		  for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
		    {
		      for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
			{	
			  for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
			    {
			      for (int kx3 = 0; kx3 < NbrSitesX; ++kx3)
				{	
				  for (int ky3 = 0; ky3 < NbrSitesY; ++ky3)
				    {
				      int kx4 = kx1 + kx2 - kx3;
				      if (kx4 < 0)
					{
					  kx4 += NbrSitesX;
					}
				      kx4 %= NbrSitesX;
				      int ky4 = ky1 + ky2 - ky3;
				      if (ky4 < 0)
					{
					  ky4 += NbrSitesY;
					}
				      ky4 %= NbrSitesY;
				      for (int i = 0; i <= 1; ++i)
					{
					  for (int j = 0; j <= 1; ++j)
					    {
					      for (int k = 0; k <= 1; ++k)
						{
						  for (int l = 0; l <= 1; ++l)
						    {						      
						      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						      CreationSigmaIndices1[NbrDensityDensityIndices] = i;
						      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						      CreationSigmaIndices2[NbrDensityDensityIndices] = j;
						      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = k;
						      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						      AnnihilationSigmaIndices2[NbrDensityDensityIndices] = l;
						      IndexLabels[NbrDensityDensityIndices] = new char[512];
						      sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d 0 %d %d %d 0 %d %d %d 0 %d %d %d 0 %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
						      ++NbrDensityDensityIndices;
						      
						      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						      CreationSigmaIndices1[NbrDensityDensityIndices] = 2 + i;
						      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						      CreationSigmaIndices2[NbrDensityDensityIndices] = 2 + j;
						      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = 2 + k;
						      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						      AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 2 + l;
						      IndexLabels[NbrDensityDensityIndices] = new char[512];
						      sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d 1 %d %d %d 1 %d %d %d 1 %d %d %d 1 %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
						      ++NbrDensityDensityIndices;
						      
						      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						      CreationSigmaIndices1[NbrDensityDensityIndices] = i;
						      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						      CreationSigmaIndices2[NbrDensityDensityIndices] = 2 + j;
						      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = k;
						      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						      AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 2 + l;
						      IndexLabels[NbrDensityDensityIndices] = new char[512];
						      sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d 0 %d %d %d 1 %d %d %d 0 %d %d %d 1 %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
						      ++NbrDensityDensityIndices;
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}	      
	    }
	  if (NbrBands == 6)
	    {
	      sprintf (FileHeader, "# kx1 ky1 spin1 sigma1 kx2 ky2 spin2 sigma2 kx3 ky3 spin3 sigma3 kx4 ky4 spin4 sigma4 <c^+ c^+ c c>");
	      NbrDensityDensityIndices = 243 * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY) * (NbrSitesX * NbrSitesY);
	      CreationMomentumIndices1 = new int[NbrDensityDensityIndices];
	      CreationSigmaIndices1 = new int[NbrDensityDensityIndices];
	      CreationMomentumIndices2 = new int[NbrDensityDensityIndices];
	      CreationSigmaIndices2 = new int[NbrDensityDensityIndices];
	      AnnihilationMomentumIndices1 = new int[NbrDensityDensityIndices];
	      AnnihilationSigmaIndices1 = new int[NbrDensityDensityIndices];
	      AnnihilationMomentumIndices2 = new int[NbrDensityDensityIndices];
	      AnnihilationSigmaIndices2 = new int[NbrDensityDensityIndices];
	      IndexLabels = new char*[NbrDensityDensityIndices];
	      NbrDensityDensityIndices = 0;
	      
	      for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
		{	
		  for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
		    {
		      for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
			{	
			  for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
			    {
			      for (int kx3 = 0; kx3 < NbrSitesX; ++kx3)
				{	
				  for (int ky3 = 0; ky3 < NbrSitesY; ++ky3)
				    {
				      int kx4 = kx1 + kx2 - kx3;
				      if (kx4 < 0)
					{
					  kx4 += NbrSitesX;
					}
				      kx4 %= NbrSitesX;
				      int ky4 = ky1 + ky2 - ky3;
				      if (ky4 < 0)
					{
					  ky4 += NbrSitesY;
					}
				      ky4 %= NbrSitesY;
				      for (int i = 0; i <= 2; ++i)
					{
					  for (int j = 0; j <= 2; ++j)
					    {
					      for (int k = 0; k <= 2; ++k)
						{
						  for (int l = 0; l <= 2; ++l)
						    {						      
						      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						      CreationSigmaIndices1[NbrDensityDensityIndices] = i;
						      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						      CreationSigmaIndices2[NbrDensityDensityIndices] = j;
						      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = k;
						      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						      AnnihilationSigmaIndices2[NbrDensityDensityIndices] = l;
						      IndexLabels[NbrDensityDensityIndices] = new char[512];
						      sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d 0 %d %d %d 0 %d %d %d 0 %d %d %d 0 %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
						      ++NbrDensityDensityIndices;
						      
						      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						      CreationSigmaIndices1[NbrDensityDensityIndices] = 3 + i;
						      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						      CreationSigmaIndices2[NbrDensityDensityIndices] = 3 + j;
						      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = 3 + k;
						      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						      AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 3 + l;
						      IndexLabels[NbrDensityDensityIndices] = new char[512];
						      sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d 1 %d %d %d 1 %d %d %d 1 %d %d %d 1 %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
						      ++NbrDensityDensityIndices;
						      
						      CreationMomentumIndices1[NbrDensityDensityIndices] = ((kx1 * NbrSitesY) + ky1);
						      CreationSigmaIndices1[NbrDensityDensityIndices] = i;
						      CreationMomentumIndices2[NbrDensityDensityIndices] = ((kx2 * NbrSitesY) + ky2);
						      CreationSigmaIndices2[NbrDensityDensityIndices] = 3 + j;
						      AnnihilationMomentumIndices1[NbrDensityDensityIndices] = ((kx3 * NbrSitesY) + ky3);
						      AnnihilationSigmaIndices1[NbrDensityDensityIndices] = k;
						      AnnihilationMomentumIndices2[NbrDensityDensityIndices] = ((kx4 * NbrSitesY) + ky4);
						      AnnihilationSigmaIndices2[NbrDensityDensityIndices] = 3 + l;
						      IndexLabels[NbrDensityDensityIndices] = new char[512];
						      sprintf(IndexLabels[NbrDensityDensityIndices], "%d %d 0 %d %d %d 1 %d %d %d 0 %d %d %d 1 %d", kx1, ky1, i, kx2, ky2, j, kx3, ky3, k, kx4, ky4, l);
						      ++NbrDensityDensityIndices;
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
      else
	{
	}

      ofstream File2;
      if (Manager.GetString("output-file") != 0)
	{
	  File2.open(Manager.GetString("output-file"), ios::binary | ios::out);
	}
      else
	{
	  char* TmpFileName = ReplaceExtensionToFileName(RightGroundStateFiles[0], "vec", "rhorho.dat");
	  if (TmpFileName == 0)
	    {
	      cout << "no vec extension was find in " << RightGroundStateFiles[0] << " file name" << endl;
	      return 0;
	    }
	  File2.open(TmpFileName, ios::binary | ios::out);
	  delete[] TmpFileName;
	}
      File2.precision(14);
      cout.precision(14);
      File2 << FileHeader << endl;
      
      if (Flag3d == false)
	{
	  if (NbrBands == 1)
	    {
	    }
	  else
	    {
	      if (NbrBands >= 2)
		{
		  for (int i = 0; i < NbrRightStates; ++i)
		    {
		      Complex TmpTotalDensity = 0.0;
		      Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		      for (int j = 0; j < NbrDensityPartialTraces; ++j)
			{
			  PartialTraces[j] = 0.0;
			}
		      int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSitesZ) + TotalKz[i];
		      for (int j = 0 ; j < NbrDensityDensityIndices; ++j)
			//		      for (int j = 0 ; j < 10; ++j)
			{
			  //			  cout << "j=" << j << ": " << CreationMomentumIndices1[j]<< " " << CreationSigmaIndices1[j]<< " " << CreationMomentumIndices2[j]<< " " << CreationSigmaIndices2[j]<< " " << AnnihilationMomentumIndices1[j]<< " " << AnnihilationSigmaIndices1[j]<< " " << AnnihilationMomentumIndices2[j]<< " " << AnnihilationSigmaIndices2[j] << endl;
			  ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], CreationMomentumIndices1[j], CreationSigmaIndices1[j], CreationMomentumIndices2[j], CreationSigmaIndices2[j], AnnihilationMomentumIndices1[j], AnnihilationSigmaIndices1[j], AnnihilationMomentumIndices2[j], AnnihilationSigmaIndices2[j]);
			  Complex TmpElement = TmpOperator.MatrixElement(RightGroundStates[i], RightGroundStates[i]);
			  File2 << IndexLabels[j] << " " << TmpElement << endl;
			}
		    }
		}
	    }
	}
      File2.close();
    }
  
  
  return 0;
}


// extract the system information from the file name
//
// stateFileName = state file name
// flag3d, flagDecoupled, flagWannier = values of flags
// nbrBands = number of bands
// nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz, totalSpin, statistics = references on system parameters
// maxBand0, maxBand1, maxBand2, maxBand3, minBand0, minBand1, minBand2, minBand3 = references on system parameters
// return value = true if no error occured

bool FTIDensityGetSystemInformationFromState(char* stateFileName, bool flag3d, bool flagDecoupled, bool flagWannier, int nbrBands,
					     int& nbrParticles, int& nbrSitesX, int& nbrSitesY, int& nbrSiteZ, int& totalKx, int& totalKy, int& totalKz, int& totalSpin, bool& statistics,
					     int& maxBand0, int& maxBand1, int& maxBand2, int& maxBand3, int& minBand0, int& minBand1, int& minBand2, int& minBand3)
{
  totalKx = 0;
  totalKy = 0;
  totalKz = 0;
  if (flag3d == false)
    {
      nbrSiteZ = 1;
      double Mass = 0.0;
      if (flagDecoupled == false)
	{
	  if(flagWannier == false) 
	    {
	      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(stateFileName,
								      nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, Mass, statistics) == false)
		{
		  return false;
		}
	      cout << stateFileName << " N=" << nbrParticles << " Nx=" << nbrSitesX << " Ny=" << nbrSitesY << " kx=" << totalKx << " ky=" << totalKy << endl;
	    }
	  else
	    {
	      if (FQHEOnSquareLatticeWannierFindSystemInfoFromVectorFileName(stateFileName,
									     nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, statistics) == false)
		{
		  return false;
		}
	    }
	}
      else
	{
	  if (nbrBands == 2)
	    {
	      if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(stateFileName,
									      nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, totalSpin, statistics) == false)
		{
		  return false;
		}
	      cout << stateFileName << " N=" << nbrParticles << " Nx=" << nbrSitesX << " Ny=" << nbrSitesY << " kx=" << totalKx << " ky=" << totalKy << " Sz=" << totalSpin << endl;
	    }
	  else
	    {
	      if (FQHEOnSquareLatticeTwoBandsWithSpinOrValleyFindSystemInfoFromVectorFileName(stateFileName,
											      nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, totalSpin, statistics) == false)
		{
		  return false;
		}
	      cout << stateFileName << " N=" << nbrParticles << " Nx=" << nbrSitesX << " Ny=" << nbrSitesY << " kx=" << totalKx << " ky=" << totalKy << " Pz=" << totalSpin << endl;
	    }
	}
      FQHEOnSquareLatticeFindMaxBandOccupationFromVectorFileName(stateFileName, maxBand0, maxBand1, maxBand2, maxBand3);
      FQHEOnSquareLatticeFindMinBandOccupationFromVectorFileName(stateFileName, minBand0, minBand1, minBand2, minBand3);
    }
  else
    {
      if (FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(stateFileName,
							     nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz, statistics) == false)
	{
	  return false;
	}
    }
  return true;
}


// get the Hilbert for any of the potential cases
//
// flag3d, flagDecoupled, flagWannier = values of flags
// nbrBands = number of bands
// nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz, totalSpin, statistics = system parameters
// maxBand0, maxBand1, maxBand2, maxBand3, minBand0, minBand1, minBand2, minBand3 = system parameters
// manager = reference on the option manager
// return value = pointer to the Hilbert space (null if an error occured)
 
ParticleOnSphere* FTIDensityGetHilbertSpace(bool flag3d, bool flagDecoupled, bool flagWannier, int nbrBands, int nbrParticles, int nbrSitesX, int nbrSitesY, int nbrSiteZ,
					    int totalKx, int totalKy, int totalKz, int totalSpin, bool statistics,
					    int maxBand0, int maxBand1, int maxBand2, int maxBand3, int minBand0, int minBand1, int minBand2, int minBand3,
					    OptionManager& manager)
{
  if (flag3d == true)
    {
      if (nbrBands == 1)
	{
	  if (statistics == true)
	    {
	      return new FermionOnCubicLatticeMomentumSpace(nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz);
	    }
	  else
	    {
	      return new BosonOnCubicLatticeMomentumSpace(nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz);
	    }
	}
      else
	{
	  if (statistics == true)
	    {
	      return new FermionOnCubicLatticeWithSpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz);
	    }
	  else
	    {
	      return new BosonOnCubicLatticeWithSU2SpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, nbrSiteZ, totalKx, totalKy, totalKz);
	    }
	 }
    }
  if (nbrBands == 1)
    {
      if (statistics == true)
	{
	  return new FermionOnSquareLatticeMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
	}
      else
	{
	  if(flagWannier == false)
	    return new BosonOnSquareLatticeMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
	  else
	    return new BosonOnSquareLatticeWannierSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKy, totalKx);
	}
      
    }
  
  if (nbrBands == 2)
    {
      if (flagDecoupled == false)
	{
	  if (statistics == true)
	    {
	      if ((maxBand0 >= 0) || (maxBand1 >= 0))
		{
		  if (maxBand0 < 0)
		    {
		      maxBand0 = 2 * nbrSitesX * nbrSitesY;
		    }
		  if (maxBand1 < 0)
		    {
		      maxBand1 = 2 * nbrSitesX * nbrSitesY;
		    }
		  if (minBand0 > (2 * nbrSitesX * nbrSitesY))
		    {
		      minBand0 = 2 * nbrSitesX * nbrSitesY;
		    }
		  if (minBand1 > (2 * nbrSitesX * nbrSitesY))
		    {
		      minBand1 = 2 * nbrSitesX * nbrSitesY;
		    }
		  if ((nbrSitesX * nbrSitesY) <= 32)
		    {
		      return new FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, minBand0, minBand1, maxBand0, maxBand1, totalKx, totalKy);
		    }
		  else
		    {
		      cout << "--max-band0 and --max-band1 options without valley and more than 32 unit cells is not yet implemented" << endl;
		      return 0;
		    }
		}
	      else
		{
		  if ((nbrSitesX * nbrSitesY) <= 32)
		    {
		      return new FermionOnSquareLatticeWithSpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
		    }
		  else
		    {
		      return new FermionOnSquareLatticeWithSpinMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
		    }
		}
	    }
	  else
	    {
	      return new BosonOnSquareLatticeWithSU2SpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
	    }
	}
      else
	{
	  if (statistics == true)
	    {
	      if ((nbrSitesX * nbrSitesY) <= 32)
		{
		  return new FermionOnSquareLatticeWithSpinMomentumSpace (nbrParticles, (totalSpin + nbrParticles) >> 1, nbrSitesX, nbrSitesY, totalKx, totalKy);
		}
	      else
		{
		  return new FermionOnSquareLatticeWithSpinMomentumSpaceLong (nbrParticles, (totalSpin + nbrParticles) >> 1, nbrSitesX, nbrSitesY, totalKx, totalKy);
		}
	    }
	  else
	    {
	      return new BosonOnSquareLatticeWithSU2SpinMomentumSpace (nbrParticles, (totalSpin + nbrParticles) >> 1, nbrSitesX, nbrSitesY, totalKx, totalKy);
	    }
 	}
    }
  if (nbrBands == 3)
    {
      if (flagDecoupled == false)
	{
	  if (statistics == true)
	    {
	      if ((maxBand0 < 0) && (maxBand1 < 0) && (maxBand2 < 0))
		{
		  if (manager.GetString("allowed-orbitals") == 0)
		    {
		      if ((nbrSitesX * nbrSitesY) <= 21)
			{
			  return new FermionOnSquareLatticeWithSU3SpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
			}
		      else
			{
			  return new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy);
			}
		    }
		  else
		    {
		      if ((nbrSitesX * nbrSitesY) <= 21)
			{
			  return new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, manager.GetString("allowed-orbitals"), totalKx, totalKy);
			}
		      else
			{
			  return new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, manager.GetString("allowed-orbitals"), totalKx, totalKy);
			}
		    }
		}
	      else
		{
		  if (maxBand0 < 0)
		    {
		      maxBand0 = nbrSitesX * nbrSitesY;
		    }
		  if (maxBand1 < 0)
		    {
		      maxBand1 = nbrSitesX * nbrSitesY;
		    }
		  if (maxBand2 < 0)
		    {
		      maxBand2 = nbrSitesX * nbrSitesY;
		    }
		  if (manager.GetString("allowed-orbitals") == 0)
		    {
		      if ((minBand0 == 0) && (minBand1 == 0) && (minBand2 == 0))
			{
			  if ((nbrSitesX * nbrSitesY) <= 21)
			    {
			      return new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, maxBand0, maxBand1, maxBand2, totalKx, totalKy);
			    }
			  else
			    {
			      return new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, maxBand0, maxBand1, maxBand2, totalKx, totalKy);
			    }
			}
		      else
			{
			  cout << minBand0<< " " << minBand1<< " " << minBand2<< " " << maxBand0<< " " << maxBand1<< " " << maxBand2 << endl;
			  if ((nbrSitesX * nbrSitesY) <= 21)
			    {
			      return new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, minBand0, minBand1, minBand2, maxBand0, maxBand1, maxBand2, totalKx, totalKy);
			    }
			  else
			    {
			      return new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, minBand0, minBand1, minBand2, maxBand0, maxBand1, maxBand2, totalKx, totalKy);
			    }
			}
		    }
		  else
		    {
		      if ((nbrSitesX * nbrSitesY) <= 21)
			{
			  return new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, manager.GetString("allowed-orbitals"), maxBand0, maxBand1, maxBand2, totalKx, totalKy);
			}
		      else
			{
			  return new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, manager.GetString("allowed-orbitals"), maxBand0, maxBand1, maxBand2, totalKx, totalKy);
			}
		    }
		}
	    }
	  else
	    {
	      cout << "bosons with 3 bands are not supported" << endl;
	      return 0;
	    }
	}
    }
  if (nbrBands == 4)
    {
      if (statistics == true)
	{
	  if ((maxBand0 < 0) && (maxBand1 < 0))
	    {
	      if ((nbrSitesX * nbrSitesY) <= 16)
		{
		  return new FermionOnSquareLatticeWithSU4SpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, totalSpin, 10000000ul);
		}
	      else
		{
		  return new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, totalSpin, 10000000ul);
		}
	    }
	  else
	    {
	      if (maxBand0 < 0)
		{
		  maxBand0 = 2 * nbrSitesX * nbrSitesY;
		}
	      if (maxBand1 < 0)
		{
		  maxBand1 = 2 * nbrSitesX * nbrSitesY;
		}
	      if ((nbrSitesX * nbrSitesY) <= 16)
		{
		  return new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, maxBand0, maxBand1, totalKx, totalKy, totalSpin, 10000000ul);
		}
	      else
		{
		  return new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, maxBand0, maxBand1, totalKx, totalKy, totalSpin, 10000000ul);
		}
	    }
	}
    }
  if (nbrBands == 6)
    {
      if (statistics == true)
	{
	  if ((nbrSitesX * nbrSitesY) <= 10)
	    {
	      return new FermionOnSquareLatticeWithSU6SpinMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, totalSpin, 10000000ul);
	    }
	  else
	    {
	      return new FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, totalKx, totalKy, totalSpin, 10000000ul);
	    }
	}
    }
  cout << "Unsupported type of Hilbert space" << endl;
  return 0;
}
 
// compute the momentum creation index with momentum transfer from the right state to the left state
//
// creationMomentumIndex = linearized
// leftTotalKx = momentum along the x direction for the left state
// leftTotalKy = momentum along the y direction for the left state
// leftTotalKz = momentum along the z direction for the left state
// totalKx = momentum along the x direction for the right state
// totalKy = momentum along the y direction for the right state
// totalKz = momentum along the z direction for the right state
// creationKx = reference on the momentum along the x direction for creation operator
// creationKy = reference on the momentum along the y direction for creation operator
// creationKz = reference on the momentum along the z direction for creation operator
// NbrSitesX = number of unit cells along the x direction
// NbrSitesY = number of unit cells along the y direction
// NbrSitesZ = number of unit cells along the z direction
// return value =momentum creation index including momentum transfer

int FTIDensityComputeMomentumTransfer(int creationMomentumIndex, int leftTotalKx, int leftTotalKy, int leftTotalKz,
				      int totalKx, int totalKy, int totalKz, int& creationKx, int& creationKy, int& creationKz,
				      int nbrSitesX, int nbrSitesY, int nbrSitesZ)
{
  int  TmpCreationMomentumIndex = creationMomentumIndex;
  creationKz = TmpCreationMomentumIndex % nbrSitesZ;
  TmpCreationMomentumIndex /= nbrSitesZ;				  
  creationKy = TmpCreationMomentumIndex % nbrSitesY;
  TmpCreationMomentumIndex /= nbrSitesY;
  creationKx = TmpCreationMomentumIndex;
  creationKx += leftTotalKx - totalKx;
  creationKy += leftTotalKy - totalKy;
  creationKz += leftTotalKz - totalKz;
  if (creationKx < 0)
    {
      creationKx += nbrSitesX;
    }
  if (creationKy < 0)
    {
      creationKy += nbrSitesY;
    }
  if (creationKz < 0)
    {
      creationKz += nbrSitesZ;
    }
  creationKx %= nbrSitesX;
  creationKy %= nbrSitesY;
  creationKz %= nbrSitesZ;
  TmpCreationMomentumIndex = (((creationKx * nbrSitesY) + creationKy) * nbrSitesZ) + creationKz;
  return TmpCreationMomentumIndex;
}
