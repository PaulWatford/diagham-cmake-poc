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
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace.h"
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
  cout.precision(14);
  
  OptionManager Manager ("FTIProjectOntoBands" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");


  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the state to be projected");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands for the input state", 1);
  (*SystemGroup) += new SingleStringOption ('\n', "allowed-orbitals", "provide an ASCII file indicating which orbitals are allowed");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name to store the projected state");
  (*OutputGroup) += new BooleanOption('\n', "normalize", "normalize the projected state");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIProjectOntoBands -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int ProjectedBandIndex = 0;
  int MaxBand0 = -1;
  int MaxBand1 = -1;
  int MaxBand2 = -1;
  int MaxBand3 = -1;
  int MinBand0 = 0;
  int MinBand1 = 0;
  int MinBand2 = 0;
  int MinBand3 = 0;
  int TotalKx = 0;
  int TotalKy = 0;
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  bool Statistics = true;
  int TotalSpin = 0;
  int NbrBands = Manager.GetInteger("nbr-subbands");

  
  if ( Manager.GetString("ground-file") == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type  FTIProjectOntoBands -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) &&  (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      char* TmpProjectedName = new char[64];
      sprintf(TmpProjectedName, "singleband_proj_band%d", ProjectedBandIndex);
      OutputFileName = ReplaceString(Manager.GetString("ground-file"), "twoband", TmpProjectedName);
      if (OutputFileName == 0)
	{
	  OutputFileName = ReplaceString(Manager.GetString("ground-file"), "threeband", TmpProjectedName);
	  if (OutputFileName == 0)
	    {
	      cout << "can't guess output name from " << Manager.GetString("ground-file") << "(should contain *band)" << endl;
	      return 0;
	    }
	}
    }
  else
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }

 if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("ground-file"), NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, Statistics) == false)
  {
     cout << "error while retrieving system parameters from file name " <<  Manager.GetString("ground-file") << endl;
     return -1;
 }
 FQHEOnSquareLatticeFindMaxBandOccupationFromVectorFileName(Manager.GetString("ground-file"), MaxBand0, MaxBand1, MaxBand2, MaxBand3);
 FQHEOnSquareLatticeFindMinBandOccupationFromVectorFileName(Manager.GetString("ground-file"), MinBand0, MinBand1, MinBand2, MinBand3);


  ParticleOnSphere* FinalSpace = 0;
  if (Statistics == true)
    {
      FinalSpace  = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
    }
  else
    {
      FinalSpace = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy); 
    }
    

  ParticleOnSphereWithSpin* InitialSpace = 0;
  if (NbrBands == 2)
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
		  InitialSpace = new FermionOnSquareLatticeWithSU2SpinAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MinBand0, MinBand1, MaxBand0, MaxBand1, TotalKx, TotalKy);
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
		  InitialSpace = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		}
	      else
		{
		  InitialSpace  = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		}
	    }
	}
      else
	{
	  InitialSpace = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy); 
	}
    }
  if (NbrBands == 3)
    {
      if (Statistics == true)
	{
	  if ((MaxBand0 < 0) && (MaxBand1 < 0) && (MaxBand2 < 0))
	    {
	      if (Manager.GetString("allowed-orbitals") == 0)
		{
		  if ((NbrSitesX * NbrSitesY) <= 21)
		    {
		      InitialSpace = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		    }
		  else
		    {
		      InitialSpace = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
		    }
		}
	      else
		{
		  if ((NbrSitesX * NbrSitesY) <= 21)
		    {
		      InitialSpace = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), TotalKx, TotalKy);
		    }
		  else
		    {
		      InitialSpace = new FermionOnSquareLatticeWithSU3SpinFilteredMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), TotalKx, TotalKy);
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
			  InitialSpace = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, MaxBand2, TotalKx, TotalKy);
			}
		      else
			{
			  InitialSpace = new FermionOnSquareLatticeWithSU3SpinAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, MaxBand0, MaxBand1, MaxBand2, TotalKx, TotalKy);
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 21)
			{
			  InitialSpace = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, MinBand0, MinBand1, MinBand2, MaxBand0, MaxBand1, MaxBand2, TotalKx, TotalKy);
			}
		      else
			{
			  InitialSpace = new FermionOnSquareLatticeWithSU3SpinAndMinMaxCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, MinBand0, MinBand1, MinBand2, MaxBand0, MaxBand1, MaxBand2, TotalKx, TotalKy);
			}
		    }
		}
	      else
		{
		  if ((NbrSitesX * NbrSitesY) <= 21)
		    {
		      InitialSpace = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), MaxBand0, MaxBand1, MaxBand2, TotalKx, TotalKy);
		    }
		  else
		    {
		      InitialSpace = new FermionOnSquareLatticeWithSU3SpinFilteredAndCapMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, Manager.GetString("allowed-orbitals"), MaxBand0, MaxBand1, MaxBand2, TotalKx, TotalKy);
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
  if (InitialSpace == 0)
    {
      cout << "unsupported number of bands" << endl;
      return 0;
    }
  
  ComplexVector InitialState;
  if (InitialState.ReadVector(Manager.GetString("ground-file")) == false)
  {
      cout << "can't open vector file " <<  Manager.GetString("ground-file") << endl;
      return -1;      
  }

  if (InitialSpace->GetLargeHilbertSpaceDimension() != InitialState.GetLargeVectorDimension())
    {
      cout << "Hilbert space dimension mismatch: " << Manager.GetString("ground-file") << " has dimension " << InitialState.GetLargeVectorDimension() << ", should be " << InitialSpace->GetLargeHilbertSpaceDimension() << endl;
      return 0;
    }
  
  ComplexVector FinalState (FinalSpace->GetHilbertSpaceDimension()); 
  InitialSpace->ProjectOntoSingleBand(&InitialState, FinalSpace, &FinalState, ProjectedBandIndex);
  cout << "weight onto band " << ProjectedBandIndex << " : " << FinalState.SqrNorm() << endl;
  if (Manager.GetBoolean("normalize") == true)
    {
      FinalState.Normalize();
    }
  FinalState.WriteVector(OutputFileName);
  

  return 0;
}
