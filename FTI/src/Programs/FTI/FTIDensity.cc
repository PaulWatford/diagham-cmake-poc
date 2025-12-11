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
  (*SystemGroup) += new BooleanOption  ('\n', "off-diagonal", "use --degenerated-groundstate as a list of states and compute all the cross density terms");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume that the FTI states are made of two decoupled FCI copies");
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

  int NbrSpaces = 1;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int* TotalKz = 0;
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  int NbrSiteZ = 0;
  int MaxBand0 = -1;
  int MaxBand1 = -1;
  int MaxBand2 = -1;
  int MaxBand3 = -1;
  int MinBand0 = 0;
  int MinBand1 = 0;
  int MinBand2 = 0;
  int MinBand3 = 0;
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
      GroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      TotalKz = new int[1];
      Coefficients = new double[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));
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
       NbrSpaces = DegeneratedFile.GetNbrLines();
       GroundStateFiles = new char* [NbrSpaces];
       TotalKx = new int[NbrSpaces];
       TotalKy = new int[NbrSpaces];
       TotalKz = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));		   
	 }
       if (DegeneratedFile.GetNbrColumns() == 1)
	 {
	   Coefficients = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] = 1.0 / ((double) NbrSpaces);
	 }
       else
	 {
	   double TmpSum = 0.0;
	   Coefficients = DegeneratedFile.GetAsDoubleArray(1);
	   for (int i = 0; i < NbrSpaces; ++i)
	     TmpSum += Coefficients[i];
	   TmpSum = 1.0 / TmpSum;
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] *= TmpSum;
	 }
    }

  if (Flag3d == false)
    {
      NbrSiteZ = 1;
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalKz[i] = 0;
	  double Mass = 0.0;
	  if (FlagDecoupled == false)
	    {
	      if(FlagWannier == false) 
		{
		  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
									  NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], Mass, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " kx=" << TotalKx[i] << " ky=" << TotalKy[i] << endl;
		}
	      else
		{
		  if (FQHEOnSquareLatticeWannierFindSystemInfoFromVectorFileName(GroundStateFiles[i],
										 NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		}
	    }
	  else
	    {
	      if (NbrBands == 2)
		{
		  if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[i],
										  NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], TotalSpin, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " kx=" << TotalKx[i] << " ky=" << TotalKy[i] << " Sz=" << TotalSpin << endl;
		}
	      else
		{
		  if (FQHEOnSquareLatticeTwoBandsWithSpinOrValleyFindSystemInfoFromVectorFileName(GroundStateFiles[i],
												  NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], TotalSpin, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " kx=" << TotalKx[i] << " ky=" << TotalKy[i] << " Pz=" << TotalSpin << endl;
		}
	    }
	  FQHEOnSquareLatticeFindMaxBandOccupationFromVectorFileName(GroundStateFiles[i], MaxBand0, MaxBand1, MaxBand2, MaxBand3);
	  FQHEOnSquareLatticeFindMinBandOccupationFromVectorFileName(GroundStateFiles[i], MinBand0, MinBand1, MinBand2, MinBand3);
	}
    }
  else
    {
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalKz[i] = 0;
	  if (FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								 NbrParticles, NbrSitesX, NbrSitesY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }

 
  GroundStates = new ComplexVector [NbrSpaces];  
  int TotalNbrSites;
  if(FlagWannier == false || (FlagWannier == true &&  TotalKx[0]>-1))
    TotalNbrSites = NbrSitesX * NbrSitesY * NbrSiteZ;
  else
    TotalNbrSites = NbrSitesY * NbrSiteZ;
  int* NbrGroundStatePerMomentumSector = new int[TotalNbrSites];
  ComplexVector** GroundStatePerMomentumSector = new ComplexVector*[TotalNbrSites];
  double** CoefficientPerMomentumSector = new double*[TotalNbrSites];
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      NbrGroundStatePerMomentumSector[i] = 0;
      GroundStatePerMomentumSector[i] = 0;
      CoefficientPerMomentumSector[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
      int TmpIndex;
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
      NbrGroundStatePerMomentumSector[TmpIndex]++; 
    }

  for (int i = 0; i < TotalNbrSites; ++i)
    {
      if (NbrGroundStatePerMomentumSector[i] > 0)
	{
	  GroundStatePerMomentumSector[i] = new ComplexVector[NbrGroundStatePerMomentumSector[i]];
	  CoefficientPerMomentumSector[i] = new double[NbrGroundStatePerMomentumSector[i]];
	}
      NbrGroundStatePerMomentumSector[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
      GroundStatePerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = GroundStates[i];
      CoefficientPerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = Coefficients[i];
      NbrGroundStatePerMomentumSector[TmpIndex]++;
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
  if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
    MaxNbrSpaces = NbrSitesX * NbrSitesY * NbrSiteZ;
  else
    MaxNbrSpaces = NbrSitesY * NbrSiteZ;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
      if (Spaces[TmpIndex] == 0)
	{
	  if (Flag3d == false)
	    {
	      if (NbrBands == 1)
		{
		  if (Statistics == true)
		    Spaces[TmpIndex] = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
		  else
		    {
		      if(FlagWannier == false)
			Spaces[TmpIndex] = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
		      else
			Spaces[TmpIndex] = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKy[i], TotalKx[i]);
		    }
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
		      if (Statistics == true)
			{
			  if ((MaxBand0 >= 0) || (MaxBand1 >= 0))
			    {
			      if (MaxBand0 < 0)
				{
				  MaxBand0 = 2 * NbrSitesX * NbrSitesY;
				}
			      if (MaxBand1 < 0)
				{
				  MaxBand1 = 2 * NbrSitesX * NbrSitesY;
				}
			      if (MinBand0 > (2 * NbrSitesX * NbrSitesY))
				{
				  MinBand0 = 2 * NbrSitesX * NbrSitesY;
				}
			      if (MinBand1 > (2 * NbrSitesX * NbrSitesY))
				{
				  MinBand1 = 2 * NbrSitesX * NbrSitesY;
				}
			      if ((NbrSitesX * NbrSitesY) <= 32)
				{
				  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MinBand0, MinBand1, MaxBand0, MaxBand1, TotalKx[i], TotalKy[i]);
				}
			      else
				{
				  cout << "--max-band0 and --max-band1 options without valley and more than 32 unit cells is not yet implemented" << endl;
				  return 0;
				}
			    }
			  else
			    {
			      if ((NbrSitesX * NbrSitesY) <= 32)
				{
				  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
				}
			      else
				{
				  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
				}
			    }
			}
		      else
			{
			  Spaces[TmpIndex] = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
			}
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
		      if (Statistics == true)
			{
			  if ((NbrSitesX * NbrSitesY) <= 32)
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
			    }
			}
		      else
			{
			  Spaces[TmpIndex] = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
			}
		      if (Manager.GetBoolean("off-diagonal") == false)
			{
			  sprintf (FileHeader, "# kx ky sigma <c^+ c>");			  
			}
		      else
			{
			  sprintf (FileHeader, "# psi_i phi_j kx ky sigma <psi_i | c^+ c | phi_j>");
			}
		      NbrDensityIndices = 2 * NbrSitesX * NbrSitesY;
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
		      if (Statistics == true)
			{
			  if ((MaxBand0 < 0) && (MaxBand1 < 0) && (MaxBand2 < 0))
			    {
			      if (Manager.GetString("allowed-orbitals") == 0)
				{
				  if ((NbrSitesX * NbrSitesY) <= 21)
				    {
				      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
				    }
				  else
				    {
				      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
				    }
				}
			      else
				{
				  if ((NbrSitesX * NbrSitesY) <= 21)
				    {
				      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), TotalKx[i], TotalKy[i]);
				    }
				  else
				    {
				      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), TotalKx[i], TotalKy[i]);
				    }
				}
			    }
			  else
			    {
			      if (MaxBand0 < 0)
				{
				  MaxBand0 = NbrSitesX * NbrSitesY;
				}
			      if (MaxBand1 < 0)
				{
				  MaxBand1 = NbrSitesX * NbrSitesY;
				}
			      if (MaxBand2 < 0)
				{
				  MaxBand2 = NbrSitesX * NbrSitesY;
				}
			      if (Manager.GetString("allowed-orbitals") == 0)
				{
				  if ((MinBand0 == 0) && (MinBand1 == 0) && (MinBand2 == 0))
				    {
				      if ((NbrSitesX * NbrSitesY) <= 21)
					{
					  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, MaxBand2, TotalKx[i], TotalKy[i]);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, MaxBand2, TotalKx[i], TotalKy[i]);
					}
				    }
				  else
				    {
				      cout << MinBand0<< " " << MinBand1<< " " << MinBand2<< " " << MaxBand0<< " " << MaxBand1<< " " << MaxBand2 << endl;
				      if ((NbrSitesX * NbrSitesY) <= 21)
					{
					  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MinBand0, MinBand1, MinBand2, MaxBand0, MaxBand1, MaxBand2, TotalKx[i], TotalKy[i]);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, MinBand0, MinBand1, MinBand2, MaxBand0, MaxBand1, MaxBand2, TotalKx[i], TotalKy[i]);
					}
				    }
				}
			      else
				{
				  if ((NbrSitesX * NbrSitesY) <= 21)
				    {
				      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), MaxBand0, MaxBand1, MaxBand2, TotalKx[i], TotalKy[i]);
				    }
				  else
				    {
				      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), MaxBand0, MaxBand1, MaxBand2, TotalKx[i], TotalKy[i]);
				    }
				}
			    }
			}
		      else
			{
			  cout << "bosons with 3 bands are not supported" << endl;
			  return 0;
			}
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
		  if (Statistics == true)
		    {
		      if ((MaxBand0 < 0) && (MaxBand1 < 0))
			{
			  if ((NbrSitesX * NbrSitesY) <= 16)
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			    }
			}
		      else
			{
			  if (MaxBand0 < 0)
			    {
			      MaxBand0 = 2 * NbrSitesX * NbrSitesY;
			    }
			  if (MaxBand1 < 0)
			    {
			      MaxBand1 = 2 * NbrSitesX * NbrSitesY;
			    }
			  if ((NbrSitesX * NbrSitesY) <= 16)
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			    }
			}
		    }
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
		  if (Statistics == true)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 10)
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU6SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			}
		      else
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			}
		    }
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
	      if (NbrBands == 1)
		{
                  if (Statistics == true)
                    Spaces[TmpIndex] = new FermionOnCubicLatticeMomentumSpace(NbrParticles, NbrSitesX, NbrSitesY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                  else
                    Spaces[TmpIndex] = new BosonOnCubicLatticeMomentumSpace(NbrParticles, NbrSitesX, NbrSitesY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                }
              else
                {
                  if (Statistics == true)
                    Spaces[TmpIndex] = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                  else
                    Spaces[TmpIndex] = new BosonOnCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                }
	    }
	}
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space (" << Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() << ") and ground state (" << GroundStates[i].GetLargeVectorDimension() << ")" << endl;
	  return 0;
	}
    }
  
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
    }
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "rho.dat");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
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
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		  Complex TmpTotalDensity = 0.0;
		  Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		  for (int j = 0; j < NbrDensityPartialTraces; ++j)
		    {
		      PartialTraces[j] = 0.0;
		    }
		  int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
		  for (int j = 0 ; j < NbrDensityIndices; ++j)
		    {
		      ParticleOnSphereDensityOperator TmpOperator ((ParticleOnSphere*) Spaces[TmpIndex], CreationMomentumIndices[j], AnnihilationMomentumIndices[j]);
		      Complex TmpElement = TmpOperator.MatrixElement(GroundStates[i], GroundStates[i]);
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
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      Complex TmpTotalDensity = 0.0;
		      Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		      for (int j = 0; j < NbrDensityPartialTraces; ++j)
			{
			  PartialTraces[j] = 0.0;
			}
		      int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
		      for (int j = 0 ; j < NbrDensityIndices; ++j)
			{
			  ParticleOnSquareLatticeWithGenericSpinBandDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], CreationMomentumIndices[j], CreationSigmaIndices[j], AnnihilationMomentumIndices[j], AnnihilationSigmaIndices[j]);
			  Complex TmpElement = TmpOperator.MatrixElement(GroundStates[i], GroundStates[i]);
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
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		}
	      cout << "Error, --off-diagonal is not implemented for a single band, no spin/valley Hilbert space" << endl;
	      return 0;
	    }
	  else
	    {
	      if (NbrBands >= 2)
		{
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      Complex TmpTotalDensity = 0.0;
		      Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		      for (int j = 0; j < NbrDensityPartialTraces; ++j)
			{
			  PartialTraces[j] = 0.0;
			}
		      int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
		      for (int j = 0; j < NbrSpaces; ++j)
			{
			  for (int k = 0 ; k < NbrDensityIndices; ++k)
			    {
			      ParticleOnSquareLatticeWithGenericSpinBandDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], CreationMomentumIndices[k], CreationSigmaIndices[k], AnnihilationMomentumIndices[k], AnnihilationSigmaIndices[k]);
			      Complex TmpElement;
			      if (Manager.GetBoolean("time-reversal") == false)
				{
				  TmpElement = TmpOperator.MatrixElement(GroundStates[i], GroundStates[j]);
				}
			      else
				{
				  TmpElement = TmpOperator.ConjugateMatrixElement(GroundStates[i], GroundStates[j]);
				}
			      if ((CreationSigmaIndices[k] == AnnihilationSigmaIndices[k]) && (CreationMomentumIndices[k] == AnnihilationMomentumIndices[k]) && (i == j))
				{
				  PartialTraces[CreationSigmaIndices[k]] += TmpElement;
				  TmpTotalDensity += TmpElement;
				}
			      File << i << " " << j << " " << IndexLabels[k] << " " << TmpElement << endl;
			    }
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
	  char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "rhorho.dat");
	  if (TmpFileName == 0)
	    {
	      cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
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
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      Complex TmpTotalDensity = 0.0;
		      Complex* PartialTraces = new Complex[NbrDensityPartialTraces];
		      for (int j = 0; j < NbrDensityPartialTraces; ++j)
			{
			  PartialTraces[j] = 0.0;
			}
		      int TmpIndex = (((TotalKx[i] * NbrSitesY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
		      for (int j = 0 ; j < NbrDensityDensityIndices; ++j)
			//		      for (int j = 0 ; j < 10; ++j)
			{
			  //			  cout << "j=" << j << ": " << CreationMomentumIndices1[j]<< " " << CreationSigmaIndices1[j]<< " " << CreationMomentumIndices2[j]<< " " << CreationSigmaIndices2[j]<< " " << AnnihilationMomentumIndices1[j]<< " " << AnnihilationSigmaIndices1[j]<< " " << AnnihilationMomentumIndices2[j]<< " " << AnnihilationSigmaIndices2[j] << endl;
			  ParticleOnSquareLatticeWithGenericSpinBandDensityDensityOperator TmpOperator ((ParticleOnSphereWithSpin*) Spaces[TmpIndex], CreationMomentumIndices1[j], CreationSigmaIndices1[j], CreationMomentumIndices2[j], CreationSigmaIndices2[j], AnnihilationMomentumIndices1[j], AnnihilationSigmaIndices1[j], AnnihilationMomentumIndices2[j], AnnihilationSigmaIndices2[j]);
			  Complex TmpElement = TmpOperator.MatrixElement(GroundStates[i], GroundStates[i]);
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
