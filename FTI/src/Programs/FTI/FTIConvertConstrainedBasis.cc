#include "Vector/ComplexVector.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"


#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

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
  OptionManager Manager ("FTIConvertConstrainedBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");


  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the state to be projected");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state (all should be defined in the same source Hilbert space)");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume that the FTI states are made of two decoupled FCI copies");
  (*SystemGroup) += new SingleStringOption ('\n', "source-allowed-orbitals", "provide an ASCII file indicating which orbitals are allowed for the source Hilbert space");
  (*SystemGroup) += new SingleStringOption ('\n', "target-allowed-orbitals", "provide an ASCII file indicating which orbitals are allowed for the target Hilbert space");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-max-band0", "maximum number of particles in band 0 for the target Hilbert space (negative if this number should be equal to the number of orbitals)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-max-band1", "maximum number of particles in band 1 for the target Hilbert space (negative if this number should be equal to the number of orbitals)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-max-band2", "maximum number of particles in band 2 for the target Hilbert space (negative if this number should be equal to the number of orbitals)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-max-band3", "maximum number of particles in band 3 for the target Hilbert space (negative if this number should be equal to the number of orbitals)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-min-band0", "minimum number of particles in band 0 for the target Hilbert space", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-min-band1", "minimum number of particles in band 1 for the target Hilbert space", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-min-band2", "minimum number of particles in band 2 for the target Hilbert space", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-min-band3", "minimum number of particles in band 3 for the target Hilbert space", 0);
  
  (*OutputGroup) += new SingleStringOption  ('\n', "output-file-name", "string that will be used in the output file names (all system parameters will be automatically added))", "projected");
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalized state once projected onto the target Hilbert space");

  (*PrecalculationGroup) += new SingleStringOption ('\n', "hilbert-directory", "directory where Hilbert spaces should be read and stored (bypassing their construction)");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIConvertConstrainedBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int TotalKx = 0;
  int TotalKy = 0;
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  int TotalSpin = 0;
  int SourceMaxBand0 = -1;
  int SourceMaxBand1 = -1;
  int SourceMaxBand2 = -1;
  int SourceMaxBand3 = -1;
  int SourceMinBand0 = -1;
  int SourceMinBand1 = -1;
  int SourceMinBand2 = -1;
  int SourceMinBand3 = -1;
  int TargetMaxBand0 = Manager.GetInteger("target-max-band0");
  int TargetMaxBand1 = Manager.GetInteger("target-max-band1");
  int TargetMaxBand2 = Manager.GetInteger("target-max-band2");
  int TargetMaxBand3 = Manager.GetInteger("target-max-band3");
  int TargetMinBand0 = Manager.GetInteger("target-min-band0");
  int TargetMinBand1 = Manager.GetInteger("target-min-band1");
  int TargetMinBand2 = Manager.GetInteger("target-min-band2");
  int TargetMinBand3 = Manager.GetInteger("target-min-band3");
  bool Statistics = true;
  int NbrBands = Manager.GetInteger("nbr-subbands");
  bool FlagDecoupled = Manager.GetBoolean("decoupled");
  int NbrSpaces = 1;
  
  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));
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
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));		   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalKx = 0;
      TotalKy = 0;
      double Mass = 0.0;
      if (FlagDecoupled == false)
	{
	  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								  NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, Mass, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " kx=" << TotalKx << " ky=" << TotalKy << endl;
	}
      else
	{
	  if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[i],
									  NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	  cout << GroundStateFiles[i] << " N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " kx=" << TotalKx << " ky=" << TotalKy << " Sz=" << TotalSpin << endl;
	}
      FQHEOnSquareLatticeFindMaxBandOccupationFromVectorFileName(GroundStateFiles[i], SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, SourceMaxBand3);
      FQHEOnSquareLatticeFindMinBandOccupationFromVectorFileName(GroundStateFiles[i], SourceMinBand0, SourceMinBand1, SourceMinBand2, SourceMinBand3);
    }
  

  char* StatisticPrefix = new char [16];
  if (Statistics == true)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  int NbrSites = NbrBands * NbrSitesX * NbrSitesY;
  char* FileSystemGeometry = new char [512];
  if (FlagDecoupled == false)
    {
      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
    }
  else
    {
      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, TotalSpin);
    }
  char* NbrBandsString = new char[16];
  switch (NbrBands)
    {
    case 1:
      sprintf (NbrBandsString, "oneband");
      break;
    case 2:
      sprintf (NbrBandsString, "twoband");
      break;
     case 3:
      sprintf (NbrBandsString, "threeband");
      break;
    case 4:
      sprintf (NbrBandsString, "fourband");
      break;
    default:
      sprintf (NbrBandsString, "zeroband");
   }
      
  char* BandMinCapPrefix = 0;
  if ((TargetMinBand0 <= 0) && (TargetMinBand1 <= 0) && (TargetMinBand2 <= 0) && (TargetMinBand3 <= 0))
    {
      if (Manager.GetString("target-allowed-orbitals") == 0)
	{
	  BandMinCapPrefix = new char [2];
	  sprintf (BandMinCapPrefix, "");
	}
      else
	{
	  BandMinCapPrefix = new char [32];
	  sprintf (BandMinCapPrefix, "_allowed_orbs");	  
	}
    }
  else
    {
      char* BandCapPrefix0 = new char [256];
      char* BandCapPrefix1 = new char [256];
      char* BandCapPrefix2 = new char [256];
      char* BandCapPrefix3 = new char [256];
      if (TargetMinBand0 <= 0) 
	{
	  sprintf (BandCapPrefix0, "");
	}
      else
	{
	  sprintf (BandCapPrefix0, "_minband0_%d", TargetMinBand0);
	}
      if (TargetMinBand1 <= 0) 
	{
	  sprintf (BandCapPrefix1, "");
	}
      else
	{
	  sprintf (BandCapPrefix1, "_minband1_%d", TargetMinBand1);
	}
      if (TargetMinBand2 <= 0) 
	{
	  sprintf (BandCapPrefix2, "");
	}
      else
	{
	  sprintf (BandCapPrefix2, "_minband2_%d", TargetMinBand2);
	}
       if (TargetMinBand3 <= 0) 
	{
	  sprintf (BandCapPrefix3, "");
	}
      else
	{
	  sprintf (BandCapPrefix3, "_minband3_%d", TargetMinBand3);
	}
     BandMinCapPrefix = new char [256 + strlen(BandCapPrefix0) + strlen(BandCapPrefix1) + strlen(BandCapPrefix2) + strlen(BandCapPrefix3)];
      if (Manager.GetString("target-allowed-orbitals") == 0)
	{
	  sprintf (BandMinCapPrefix, "%s%s%s%s", BandCapPrefix0, BandCapPrefix1, BandCapPrefix2, BandCapPrefix3);
	}
      else
	{
	  sprintf (BandMinCapPrefix, "_allowed_orbs%s%s%s%s", BandCapPrefix0, BandCapPrefix1, BandCapPrefix2, BandCapPrefix3);
	}
    }
  char* BandMaxCapPrefix = 0;
  if ((TargetMaxBand0 < 0) && (TargetMaxBand1 < 0) && (TargetMaxBand2 < 0))
    {
      if (Manager.GetString("target-allowed-orbitals") == 0)
	{
	  BandMaxCapPrefix = new char [2];
	  sprintf (BandMaxCapPrefix, "");
	}
      else
	{
	  BandMaxCapPrefix = new char [32];
	  sprintf (BandMaxCapPrefix, "_allowed_orbs");	  
	}
    }
  else
    {
      char* BandCapPrefix0 = new char [256];
      char* BandCapPrefix1 = new char [256];
      char* BandCapPrefix2 = new char [256];
      char* BandCapPrefix3 = new char [256];
      if (TargetMaxBand0 < 0) 
	{
	  sprintf (BandCapPrefix0, "");
	}
      else
	{
	  sprintf (BandCapPrefix0, "_maxband0_%d", TargetMaxBand0);
	}
      if (TargetMaxBand1 < 0) 
	{
	  sprintf (BandCapPrefix1, "");
	}
      else
	{
	  sprintf (BandCapPrefix1, "_maxband1_%d", TargetMaxBand1);
	}
      if (TargetMaxBand2 < 0) 
	{
	  sprintf (BandCapPrefix2, "");
	}
      else
	{
	  sprintf (BandCapPrefix2, "_maxband2_%d", TargetMaxBand2);
	}
       if (TargetMaxBand3 < 0) 
	{
	  sprintf (BandCapPrefix3, "");
	}
      else
	{
	  sprintf (BandCapPrefix3, "_maxband3_%d", TargetMaxBand3);
	}
     BandMaxCapPrefix = new char [256 + strlen(BandCapPrefix0) + strlen(BandCapPrefix1) + strlen(BandCapPrefix2) + strlen(BandCapPrefix3)];
      if ((Manager.GetString("target-allowed-orbitals") == 0) || ((TargetMinBand0 <= 0) && (TargetMinBand1 <= 0) && (TargetMinBand2 <= 0) && (TargetMinBand3 <= 0)))
	{
	  sprintf (BandMaxCapPrefix, "%s%s%s%s", BandCapPrefix0, BandCapPrefix1, BandCapPrefix2, BandCapPrefix3);
	}
      else
	{
	  sprintf (BandMaxCapPrefix, "_allowed_orbs%s%s%s%s", BandCapPrefix0, BandCapPrefix1, BandCapPrefix2, BandCapPrefix3);
	}
    }
  char* FilePrefix = new char [512 + strlen(FileSystemGeometry) + strlen(BandMaxCapPrefix) + strlen(BandMinCapPrefix) +  strlen(Manager.GetString("output-file-name"))];
  sprintf (FilePrefix, "%s_%s_%s%s%s_%s", StatisticPrefix, NbrBandsString, Manager.GetString("output-file-name"), BandMinCapPrefix, BandMaxCapPrefix, FileSystemGeometry);

  
  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type  FTIConvertConstrainedBasis -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) &&  (IsFile(Manager.GetString("ground-file")) == false))
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
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }


  ComplexVector* InitialStates = new ComplexVector [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (InitialStates[i].ReadVector(GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i]  << endl;
	  return 0;      
	}
    }

  ParticleOnSphereWithSpin* InputSpace = 0;
  ParticleOnSphereWithSpin* OutputSpace = 0;
  
  if (NbrBands == 1)
    {
      cout << "one band is not supported" << endl;
      return 0;
    }
  if (NbrBands == 2)
    {
      if (FlagDecoupled == false)
	{
	  if (Statistics == true)
	    {
	      if ((NbrSitesX * NbrSitesY) <= 32)
		{
		  InputSpace = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		}
	      else
		{
		  InputSpace = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		}
	    }
	  else
	    {
	      InputSpace = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
	    }
	}
      else
	{
	  if (Statistics == true)
	    {
	      if ((NbrSitesX * NbrSitesY) <= 32)
		{
		  InputSpace = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		}
	      else
		{
		  InputSpace = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		}
	    }
	  else
	    {
	      InputSpace = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
	    }
	}
    }
  if (NbrBands == 3)
    {
      if (FlagDecoupled == false)
	{
	  if (Statistics == true)
	    {
	      if ((SourceMaxBand0 < 0) && (SourceMaxBand1 < 0) && (SourceMaxBand2 < 0))
		{
		  if (Manager.GetString("source-allowed-orbitals") == 0)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  InputSpace = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
			}
		      else
			{
			  InputSpace = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  InputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("source-allowed-orbitals"), TotalKx, TotalKy);
			}
		      else
			{
			  InputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("source-allowed-orbitals"), TotalKx, TotalKy);
			}
		    }
		}
	      else
		{
		  if (SourceMaxBand0 < 0)
		    {
		      SourceMaxBand0 = NbrSitesX * NbrSitesY;
		    }
		  if (SourceMaxBand1 < 0)
		    {
		      SourceMaxBand1 = NbrSitesX * NbrSitesY;
		    }
		  if (SourceMaxBand2 < 0)
		    {
		      SourceMaxBand2 = NbrSitesX * NbrSitesY;
		    }
		  if (Manager.GetString("source-allowed-orbitals") == 0)
		    {
		      if ((SourceMinBand0 == 0) && (SourceMinBand1 == 0) && (SourceMinBand2 == 0))
			{
			  if ((NbrSitesX * NbrSitesY) <= 21)
			    {
			      InputSpace = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, TotalKx, TotalKy);
			    }
			  else
			    {
			      InputSpace = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, TotalKx, TotalKy);
			    }
			}
		      else
			{
			  if ((NbrSitesX * NbrSitesY) <= 21)
			    {
			      InputSpace = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, SourceMinBand0, SourceMinBand1, SourceMinBand2, SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, TotalKx, TotalKy);
			    }
			  else
			    {
			      InputSpace = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, SourceMinBand0, SourceMinBand1, SourceMinBand2, SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, TotalKx, TotalKy);
			    }
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  InputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("source-allowed-orbitals"), SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, TotalKx, TotalKy);
			}
		      else
			{
			  InputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("source-allowed-orbitals"), SourceMaxBand0, SourceMaxBand1, SourceMaxBand2, TotalKx, TotalKy);
			}
		    }
		}

	      if ((TargetMaxBand0 < 0) && (TargetMaxBand1 < 0) && (TargetMaxBand2 < 0))
		{
		  if (Manager.GetString("target-allowed-orbitals") == 0)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  OutputSpace = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
			}
		      else
			{
			  OutputSpace = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  OutputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("target-allowed-orbitals"), TotalKx, TotalKy);
			}
		      else
			{
			  OutputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("target-allowed-orbitals"), TotalKx, TotalKy);
			}
		    }
		}
	      else
		{
		  if (TargetMaxBand0 < 0)
		    {
		      TargetMaxBand0 = NbrSitesX * NbrSitesY;
		    }
		  if (TargetMaxBand1 < 0)
		    {
		      TargetMaxBand1 = NbrSitesX * NbrSitesY;
		    }
		  if (TargetMaxBand2 < 0)
		    {
		      TargetMaxBand2 = NbrSitesX * NbrSitesY;
		    }
		  if (Manager.GetString("target-allowed-orbitals") == 0)
		    {
		      if ((TargetMinBand0 == 0) && (TargetMinBand1 == 0) && (TargetMinBand2 == 0))
			{
			  if ((NbrSitesX * NbrSitesY) <= 21)
			    {
			      OutputSpace = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TargetMaxBand0, TargetMaxBand1, TargetMaxBand2, TotalKx, TotalKy);
			    }
			  else
			    {
			      OutputSpace = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TargetMaxBand0, TargetMaxBand1, TargetMaxBand2, TotalKx, TotalKy);
			    }
			}
		      else
			{
			  if ((NbrSitesX * NbrSitesY) <= 21)
			    {
			      OutputSpace = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TargetMinBand0, TargetMinBand1, TargetMinBand2, TargetMaxBand0, TargetMaxBand1, TargetMaxBand2, TotalKx, TotalKy);
			    }
			  else
			    {
			      OutputSpace = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TargetMinBand0, TargetMinBand1, TargetMinBand2, TargetMaxBand0, TargetMaxBand1, TargetMaxBand2, TotalKx, TotalKy);
			    }
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  OutputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("target-allowed-orbitals"), TargetMaxBand0, TargetMaxBand1, TargetMaxBand2, TotalKx, TotalKy);
			}
		      else
			{
			  OutputSpace = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("target-allowed-orbitals"), TargetMaxBand0, TargetMaxBand1, TargetMaxBand2, TotalKx, TotalKy);
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
  if (NbrBands == 4)
    {
      if (Statistics == true)
	{
	  if ((SourceMaxBand0 < 0) && (SourceMaxBand1 < 0))
	    {
	      if ((NbrSitesX * NbrSitesY) <= 16)
		{
		  InputSpace = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, 10000000ul);
		}
	      else
		{
		  InputSpace = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, 10000000ul);
		}
	    }
	  else
	    {
	      if (SourceMaxBand0 < 0)
		{
		  SourceMaxBand0 = 2 * NbrSitesX * NbrSitesY;
		}
	      if (SourceMaxBand1 < 0)
		{
		  SourceMaxBand1 = 2 * NbrSitesX * NbrSitesY;
		}
	      if ((NbrSitesX * NbrSitesY) <= 16)
		{
		  InputSpace = new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, SourceMaxBand0, SourceMaxBand1, TotalKx, TotalKy, TotalSpin, Manager.GetString("hilbert-directory"), 10000000ul);
		}
	      else
		{
		  InputSpace = new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, SourceMaxBand0, SourceMaxBand1, TotalKx, TotalKy, TotalSpin, Manager.GetString("hilbert-directory"), 10000000ul);
		}
	    }

	  if ((TargetMaxBand0 < 0) && (TargetMaxBand1 < 0))
	    {
	      if ((NbrSitesX * NbrSitesY) <= 16)
		{
		  OutputSpace = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, 10000000ul);
		}
	      else
		{
		  OutputSpace = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, 10000000ul);
		}
	    }
	  else
	    {
	      if (TargetMaxBand0 < 0)
		{
		  TargetMaxBand0 = 2 * NbrSitesX * NbrSitesY;
		}
	      if (TargetMaxBand1 < 0)
		{
		  TargetMaxBand1 = 2 * NbrSitesX * NbrSitesY;
		}
	      if ((NbrSitesX * NbrSitesY) <= 16)
		{
		  OutputSpace = new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TargetMaxBand0, TargetMaxBand1, TotalKx, TotalKy, TotalSpin, Manager.GetString("hilbert-directory"), 10000000ul);
		}
	      else
		{
		  OutputSpace = new FermionOnSquareLatticeWithSU4SpinAndValleyCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TargetMaxBand0, TargetMaxBand1, TotalKx, TotalKy, TotalSpin, Manager.GetString("hilbert-directory"), 10000000ul);
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
	      InputSpace = new FermionOnSquareLatticeWithSU6SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, 10000000ul);
	    }
	  else
	    {
	      InputSpace = new FermionOnSquareLatticeWithSU6SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, 10000000ul);
	    }
	}
    }

  ComplexVector TmpOutputState;
  char* TmpExtention = new char[256];
  if (FlagDecoupled == false)
    {
      sprintf (TmpExtention, "_kx_%d_ky_%d", TotalKx, TotalKy);
    }
  else
    {
      sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", TotalKx, TotalKy, TotalSpin);
    }
 for (int i = 0; i < NbrSpaces; ++i)
    {
      if (InitialStates[i].GetVectorDimension() == InputSpace->GetHilbertSpaceDimension())
	{
	  TmpOutputState = OutputSpace->ConvertFromNbodyBasis(InitialStates[i], *InputSpace);
	  double TmpNorm = TmpOutputState.Norm();
	  cout.precision(14);
	  cout << "Weight of projected " << GroundStateFiles[i] << " : " << (TmpNorm * TmpNorm) << endl;
	  if (Manager.GetBoolean("normalize"))
	    {
	      if (TmpNorm != 0.0)
		{
		  TmpOutputState /= TmpNorm;
		}
	      else
		{
		  cout << "warning, projected " << GroundStateFiles[i] << " state is a null vector" << endl;
		}
	    }
	  char* TmpOutputFile = new char[strlen(FilePrefix) + strlen(TmpExtention) + 16];
	  sprintf (TmpOutputFile, "%s%s.%d.vec", FilePrefix, TmpExtention, i);
	  if (TmpOutputState.WriteVector(TmpOutputFile) == false)
	    {
	      cout << "can't write vector file " << TmpOutputFile  << endl;
	      return 0;      
	    }
	  delete[] TmpOutputFile;
	}
      else
	{
	  cout << "dimension mismatch for state " << GroundStateFiles[i] << "(is " << InitialStates[i].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
	}
    }

  return 0;
}
