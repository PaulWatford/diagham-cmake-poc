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

#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU6SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU12SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU12SpinMomentumSpaceLong.h"

#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"

#include "Operator/ParticleOnSquareLatticeWithGenericSpinBandDensityOperator.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume that the FTI states are made of two decoupled FCI copies");
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption  ('\n', "Wannier", "Wannier basis");
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
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int NbrSiteZ = 0;
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
									  NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Mass, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSiteX << " Ny=" << NbrSiteY << " kx=" << TotalKx[i] << " ky=" << TotalKy[i] << endl;
		}
	      else
		{
		  if (FQHEOnSquareLatticeWannierFindSystemInfoFromVectorFileName(GroundStateFiles[i],
										 NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Statistics) == false)
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
										  NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSiteX << " Ny=" << NbrSiteY << " kx=" << TotalKx[i] << " ky=" << TotalKy[i] << " Sz=" << TotalSpin << endl;
		}
	      else
		{
		  if (FQHEOnSquareLatticeTwoBandsWithSpinOrValleyFindSystemInfoFromVectorFileName(GroundStateFiles[i],
												  NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSiteX << " Ny=" << NbrSiteY << " kx=" << TotalKx[i] << " ky=" << TotalKy[i] << " Pz=" << TotalSpin << endl;
		}
	    }
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
								 NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }

 
  GroundStates = new ComplexVector [NbrSpaces];  
  int TotalNbrSites;
  if(FlagWannier == false || (FlagWannier == true &&  TotalKx[0]>-1))
    TotalNbrSites = NbrSiteX * NbrSiteY * NbrSiteZ;
  else
    TotalNbrSites = NbrSiteY * NbrSiteZ;
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
	TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
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
	TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
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
    MaxNbrSpaces = NbrSiteX * NbrSiteY * NbrSiteZ;
  else
    MaxNbrSpaces = NbrSiteY * NbrSiteZ;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
      if (Spaces[TmpIndex] == 0)
	{
	  if (Flag3d == false)
	    {
	      if (NbrBands == 1)
		{
		  if (Statistics == true)
		    Spaces[TmpIndex] = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		  else
		    {
		      if(FlagWannier == false)
			Spaces[TmpIndex] = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		      else
			Spaces[TmpIndex] = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKy[i], TotalKx[i]);
		    }
		}
	      if (NbrBands == 2)
		{
		  if (FlagDecoupled == false)
		    {
		      if (Statistics == true)
			{
			  if ((NbrSiteX * NbrSiteY) <= 32)
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			    }
			}
		      else
			{
			  Spaces[TmpIndex] = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			}
		      sprintf (FileHeader, "# kx ky sigma sigma' <c^+ c>");
		      NbrDensityIndices = 4 * NbrSiteX * NbrSiteY;
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
		      for (int kx = 0; kx < NbrSiteX; ++kx)
			{	
			  for (int ky = 0; ky < NbrSiteY; ++ky)
			    {
			      for (int i = 0; i <= 1; ++i)
				{
				  for (int j = 0; j <= 1; ++j)
				    {
				      CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				      AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
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
			  if ((NbrSiteX * NbrSiteY) <= 32)
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			    }
			}
		      else
			{
			  Spaces[TmpIndex] = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			}
		      sprintf (FileHeader, "# kx ky sigma <c^+ c>");
		      NbrDensityIndices = 2 * NbrSiteX * NbrSiteY;
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
		      for (int kx = 0; kx < NbrSiteX; ++kx)
			{	
			  for (int ky = 0; ky < NbrSiteY; ++ky)
			    {
			      for (int i = 0; i <= 1; ++i)
				{
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
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
			  if ((NbrSiteX * NbrSiteY) <= 32)
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
			    }
			}
		      else
			{
			  cout << "bosons with 3 bands are not supported" << endl;
			  return 0;
			}
		      sprintf (FileHeader, "# kx ky sigma sigma' <c^+ c>");
		      NbrDensityIndices = 9 * NbrSiteX * NbrSiteY;
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
		      for (int kx = 0; kx < NbrSiteX; ++kx)
			{	
			  for (int ky = 0; ky < NbrSiteY; ++ky)
			    {
			      for (int i = 0; i <= 2; ++i)
				{
				  for (int j = 0; j <= 2; ++j)
				    {
				      CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				      AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
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
		      if ((NbrSiteX * NbrSiteY) <= 16)
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			}
		      else
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			}
		    }
		  sprintf (FileHeader, "# kx ky spin sigma sigma' <c^+ c>");
		  NbrDensityIndices = 8 * NbrSiteX * NbrSiteY;
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
		  for (int kx = 0; kx < NbrSiteX; ++kx)
		    {	
		      for (int ky = 0; ky < NbrSiteY; ++ky)
			{
			  for (int i = 0; i <= 1; ++i)
			    {
			      for (int j = 0; j <= 1; ++j)
				{
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = j;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d 0 %d %d", kx, ky, i, j);
				  ++NbrDensityIndices;
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
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
		      if ((NbrSiteX * NbrSiteY) <= 16)
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU6SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			}
		      else
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, 10000000ul);
			}
		    }
		  sprintf (FileHeader, "# kx ky spin sigma sigma' <c^+ c>");
		  NbrDensityIndices = 18 * NbrSiteX * NbrSiteY;
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
		  for (int kx = 0; kx < NbrSiteX; ++kx)
		    {	
		      for (int ky = 0; ky < NbrSiteY; ++ky)
			{
			  for (int i = 0; i <= 2; ++i)
			    {
			      for (int j = 0; j <= 2; ++j)
				{
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  CreationSigmaIndices[NbrDensityIndices] = i;
				  AnnihilationSigmaIndices[NbrDensityIndices] = j;
				  IndexLabels[NbrDensityIndices] = new char[256];
				  sprintf(IndexLabels[NbrDensityIndices], "%d %d 0 %d %d", kx, ky, i, j);
				  ++NbrDensityIndices;
				  CreationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
				  AnnihilationMomentumIndices[NbrDensityIndices] = ((kx * NbrSiteY) + ky);
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
                    Spaces[TmpIndex] = new FermionOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                  else
                    Spaces[TmpIndex] = new BosonOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                }
              else
                {
                  if (Statistics == true)
                    Spaces[TmpIndex] = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                  else
                    Spaces[TmpIndex] = new BosonOnCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
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
		  int TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
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
  File.close();

  return 0;
}
