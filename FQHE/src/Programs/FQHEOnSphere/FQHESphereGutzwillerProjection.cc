#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzGutzwillerProjection.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereGutzwillerProjection" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "name of the file on which the Gutzwiller projection has to be applied");
  (*SystemGroup) += new SingleStringOption  ('\n', "file-list", "single column file describing a list to states that have to be projected (should all have the same quantum numbers)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (while appending .x.vec at the end of each stored vector)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereGutzwillerProjection -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
  int NbrSpaces = 1;
  ComplexVector* InputStates = 0;
  char** InputStateFiles = 0;
  int TotalLz = 0;
  int TotalSz = 0;
  int LzSymmetry = 0;
  int SzSymmetry = 0;
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool LzSymmetrizedBasis =false ;
  bool LzMinusParity = false;
  int NbrParticles = 0;
  int LzMax = 0;
  bool FermionFlag = true;
  bool AllSzFlag = false;
  bool SU2SpinFlag = true;
  bool SU3SpinFlag = false;
  bool SU4SpinFlag = false;
  bool GutzwillerFlag = false;
  unsigned long MemorySpace = 9ul << 20;
  
  if ((Manager.GetString("input-file") == 0) && (Manager.GetString("file-list") == 0))
    {
      cout << "error, an input state file should be provided. See man page for option syntax or type FQHESphereGutzwillerProjection -h" << endl;
      return -1;
    }
  if ((Manager.GetString("input-file") != 0) && 
      (IsFile(Manager.GetString("input-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
      return -1;
    }
  if ((Manager.GetString("file-list") != 0) && 
      (IsFile(Manager.GetString("file-list")) == false))
    {
      cout << "can't open file " << Manager.GetString("file-list") << endl;
      return -1;
    }

  if (Manager.GetString("file-list") == 0)
    {
      InputStateFiles = new char* [1];
      InputStateFiles[0] = new char [strlen(Manager.GetString("input-file")) + 1];
      strcpy (InputStateFiles[0], Manager.GetString("input-file"));
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("file-list")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       InputStateFiles = new char* [NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   InputStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (InputStateFiles[i], DegeneratedFile(0, i));		   
	 }
    }

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(InputStateFiles[0], NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetry, LzSymmetry, FermionFlag, AllSzFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " << InputStateFiles[0] << endl;
      return -1;
      
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


  ParticleOnSphereWithSpin* Space = 0;
  ParticleOnSphereWithSpin* ProjectedSpace = 0;

  if (FermionFlag == true)
    {
      if (AllSzFlag)
	{
	  if (LzSymmetrizedBasis == false)
	    {
	      if (SzSymmetrizedBasis == false)
		{
		  Space = new FermionOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, MemorySpace);
		  ProjectedSpace = new FermionOnSphereWithSpinAllSzGutzwillerProjection(NbrParticles, TotalLz, LzMax, MemorySpace);
		}
	      else
		{
		  Space = new FermionOnSphereWithSpinAllSzSzSymmetry(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
		  ProjectedSpace = new FermionOnSphereWithSpinAllSzGutzwillerProjectionSzSymmetry (NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
		}
	    }
	  else
	    {
	      cout << "Gutzwiller projection is not implemented for systems with inversion symmetry Lz<->-Lz" << endl;
	      return 0;
	      // if (SzSymmetrizedBasis == false)
	      // 	{
	      // 	  Space = new FermionOnSphereWithSpinAllSzLzSymmetry (NbrParticles, LzMax, LzMinusParity, MemorySpace);
	      // 	}
	      // else
	      // 	{
	      // 	  Space =  new FermionOnSphereWithSpinAllSzLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, MemorySpace);
	      // 	}
	    }
	}
      else
	{
	  cout << "Gutzwiller projection is not implemented for systems with a fixed total Sz" << endl;
	  return 0;
// 	  if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
// 	    {
// #ifdef __64_BITS__
// 	      if (LzMax <= 31)
// #else
// 		if (LzMax <= 15)
// #endif
// 		  {
// 		    Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
// 		  }
// 		else
// 		  {
// #ifdef __128_BIT_LONGLONG__
// 		    if (LzMax <= 63)
// #else
// 		      if (LzMax <= 31)
// #endif
// 			{
// 			  Space = new FermionOnSphereWithSpinLong(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
// 			}
// 		      else
// 			{
// 			  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
// 			  return -1;
// 			}	
// 		  }
// 	    }
// 	  else
// 	    {
// #ifdef __128_BIT_LONGLONG__
// 	      if (LzMax >= 61)
// #else
// 		if (LzMax >= 29)
// #endif
// 		  {
// 		    cout << "States of this Hilbert space cannot be represented in a single word." << endl;
// 		    return -1;
// 		  }	
// 	      if (SzSymmetrizedBasis == true) 
// 		if (LzSymmetrizedBasis == false)
// 		  {
// #ifdef __64_BITS__
// 		    if (LzMax <= 28)
// #else
// 		      if (LzMax <= 13)
// #endif
// 			{
// 			  if (Manager.GetString("load-hilbert") == 0)
// 			    Space = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
// 			  else
// 			    Space = new FermionOnSphereWithSpinSzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
// 			}
// 		      else
// 			{
// 			  if (Manager.GetString("load-hilbert") == 0)
// 			    Space = new FermionOnSphereWithSpinSzSymmetryLong(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
// 			  else
// 			    Space = new FermionOnSphereWithSpinSzSymmetryLong(Manager.GetString("load-hilbert"), MemorySpace);
// 			}
// 		  }
// 		else
// #ifdef __64_BITS__
// 		  if (LzMax <= 28)
// #else
// 		    if (LzMax <= 13)
// #endif
// 		      {
// 			if (Manager.GetString("load-hilbert") == 0)
// 			  {
// 			    Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, LzMax, SzMinusParity,
// 									    LzMinusParity, MemorySpace);
// 			  }
// 			else
// 			  Space = new FermionOnSphereWithSpinLzSzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
// 		      }
// 		    else
// 		      {
// 			if (Manager.GetString("load-hilbert") == 0)
// 			  {
// 			    Space = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrParticles, LzMax, SzMinusParity,
// 										LzMinusParity, MemorySpace);
// 			  }
// 			else
// 			  Space = new FermionOnSphereWithSpinLzSzSymmetryLong(Manager.GetString("load-hilbert"), MemorySpace);
			
// 		      }
// 		  else
// #ifdef __64_BITS__
// 		    if (LzMax <= 28)
// #else
// 		      if (LzMax <= 13)
// #endif
// 			{
// 			  if (Manager.GetString("load-hilbert") == 0)
// 			    Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
// 			  else
// 			    Space = new FermionOnSphereWithSpinLzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);	      
// 			}
// 		      else
// 			{
// 			  if (Manager.GetString("load-hilbert") == 0)
// 			    Space = new FermionOnSphereWithSpinLzSymmetryLong(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
// 			  else
// 			    Space = new FermionOnSphereWithSpinLzSymmetryLong(Manager.GetString("load-hilbert"), MemorySpace);	      
// 			}
// 	    }
// 	}
	}
    }
  else
    {
      cout << "gutzwiller projection is not implemented for bosons" << endl;
      return 0;
    }
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      cout << "projecting " << InputStateFiles[i] << endl;
      RealVector TmpVector;
      if (TmpVector.ReadVector (InputStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << InputStateFiles[i] << endl;
	  return -1;      
	}
      if (TmpVector.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << InputStateFiles[i] << " has the wrong dimension (is " << TmpVector.GetVectorDimension() << ", should be " << Space->GetHilbertSpaceDimension() << ")" << endl;
	  return 0;
	}
      RealVector TmpTargetVector = ProjectedSpace->GutzwillerProjection(TmpVector, Space);

      
      double TmpWeight = TmpTargetVector.SqrNorm();
      cout << "   weight of the projected state = " << TmpWeight << endl;
      if (TmpWeight > MACHINE_PRECISION)
	{
	  TmpTargetVector /= sqrt(TmpWeight);
	  char* TmpOutputName;
	  if (Manager.GetString("output-file") != 0)
	    {
	      TmpOutputName = new char[strlen(Manager.GetString("output-file"))+ 24];
	      sprintf(TmpOutputName, "%s.%d.vec", Manager.GetString("output-file"), i);
	    }
	  else
	    {
	      char* TmpString = 0;
	      TmpOutputName = ReplaceString(InputStateFiles[i], "_su2_", "_gutzwiller_projected_su2_");
	      if (TmpOutputName == 0)
		{
		  cout << "cannot build output file name from file name " << InputStateFiles[i] << endl;
		  return 0;
		}
	    }
	  if (TmpTargetVector.WriteVector(TmpOutputName) == false)
	    {
	      cout << "can't write " << TmpOutputName << endl;
	      return 0;
	    }
	}
    }

  delete Space;
  delete ProjectedSpace;
  return 0;
}
