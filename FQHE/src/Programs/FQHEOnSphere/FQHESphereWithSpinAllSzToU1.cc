#include "config.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsAllSzToU1" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file that contains the SU2WithTunneling state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('t', "tunneling-amp", "tunneling amplitude", 0.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsAllSzToU1 -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "no input state " << endl << "see man page for option syntax or type FQHESphereFermionsAllSzToU1 -h" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  //int PairParity = Manager.GetInteger("pair-parity");

  double tunneling = Manager.GetDouble("tunneling-amp");

  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  bool FermionFlag = false;

//  if (NbrParticles == 0)
//  if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["state"])->GetString(), NbrParticles, LzMax, TotalLz, 
//                    TotalSz, FermionFlag) == false)
//    {
//      cout << "error while retrieving system informations from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
//      return -1;
//    }

  char* StateFileName = ((SingleStringOption*) Manager["state"])->GetString();

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
            cout << "Error " << Manager.GetString("statistics") << " is an undefined statistics" << endl;
          }
       }
    }
  else
   FermionFlag = true;  

  if (FermionFlag)
    cout << Manager.GetString("statistics") << " N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << endl;
  else
    cout << Manager.GetString("statistics") << " N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << endl;

  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  RealVector State;
  if (State.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }


  unsigned long MemorySpace = 9l << 20;
  char* OutputName = new char [512 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  if (FermionFlag)
     sprintf (OutputName, "fermions_%s_n_%d_2s_%d_t_%f_lz_%d.0.vec", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, tunneling, TotalLz);
  else
    sprintf (OutputName, "bosons_%s_n_%d_2s_%d_t_%f_lz_%d.0.vec", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, tunneling, TotalLz); 

  if (FermionFlag == true)
    {
      FermionOnSphere* U1Space = 0;
#ifdef __64_BITS__
      if (LzMax <= 63)
#else
	if (LzMax <= 31)
#endif
	  {
	    U1Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
	  }
	else
	  {
	    cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	    return -1;
	  }	

	  FermionOnSphereWithSpinAllSz* Space;
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
    
	  RealVector OutputState = Space->ForgeU1FromTunneling(State, *U1Space);
	  OutputState.WriteVector(OutputName);
	
	delete Space;
      delete U1Space;
    }
 else //.....bosons.....
   {
    
      if (LzSymmetrizedBasis == false)
       {
         BosonOnSphere* U1Space = new BosonOnSphere(NbrParticles, TotalLz, LzMax);

         BosonOnSphereWithSpinAllSz* Space = new BosonOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);

         //int PairParity = -1;
         //if ( PairParity >=0 ) 
         //   Space = new BosonOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, PairParity, MemorySpace);
         // else
         //   Space = new BosonOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);
    
         RealVector OutputState = Space->ForgeU1FromTunneling(State, *U1Space);
         OutputState.WriteVector(OutputName);  
         delete Space;
         delete U1Space;
       }
     else
      {
        cout << "Lz-symmetrized states not available for Bosons with Spin."<<endl;
        return -1;
      }
     
   }   
  return 0;
}

