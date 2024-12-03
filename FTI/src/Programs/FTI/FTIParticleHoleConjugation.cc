#include "Vector/ComplexVector.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"


#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"


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
  
  OptionManager Manager ("FTIParticleHoleConjugation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");


  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the state to be projected");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name to store the projected state");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIParticleHoleConjugation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int TotalKx = 0;
  int TotalKy = 0;
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  bool Statistics = true;
  int TotalSpin = 0;
  
  if ( Manager.GetString("ground-file") == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type  FTIParticleHoleConjugation -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) &&  (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }

  if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("ground-file"), NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, TotalSpin, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " <<  Manager.GetString("ground-file") << endl;
      return -1;
    }

  int TargetTotalKx = (NbrSitesX - TotalKx) % NbrSitesX;
  int TargetTotalKy = (NbrSitesY - TotalKy) % NbrSitesY;
  int TargetNbrParticles = (NbrSitesX * NbrSitesY) - NbrParticles;
  
  if (Statistics == false)
    {
      cout << "particle-hole conjugation only works for fermionic states" << endl;
      return 0;
    }
  
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      char* TmpProjectedName1 = new char[64];
      char* TmpProjectedName2 = new char[64];
      char* TmpName = ReplaceString(Manager.GetString("ground-file"), "fermions", "fermions_ph");
      sprintf(TmpProjectedName1, "n_%d_", NbrParticles);
      sprintf(TmpProjectedName2, "n_%d_", TargetNbrParticles);
      char* TmpName2 = ReplaceString(TmpName, TmpProjectedName1, TmpProjectedName2);
      sprintf(TmpProjectedName1, "kx_%d",TotalKx);
      sprintf(TmpProjectedName2, "kx_%d", TargetTotalKx);
      char* TmpName3 = ReplaceString(TmpName2, TmpProjectedName1, TmpProjectedName2);
      sprintf(TmpProjectedName1, "ky_%d",TotalKy);
      sprintf(TmpProjectedName2, "ky_%d", TargetTotalKy);
      OutputFileName = ReplaceString(TmpName3, TmpProjectedName1, TmpProjectedName2);
    }
  else
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }

  FermionOnSquareLatticeMomentumSpace* TargetSpace = new FermionOnSquareLatticeMomentumSpace (TargetNbrParticles, NbrSitesX, NbrSitesY, TargetTotalKx, TargetTotalKy);
  FermionOnSquareLatticeMomentumSpace* InitialSpace = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
  
  ComplexVector InitialState;
  if (InitialState.ReadVector(Manager.GetString("ground-file")) == false)
  {
      cout << "can't open vector file " <<  Manager.GetString("ground-file") << endl;
      return -1;      
  }
  if (InitialState.GetVectorDimension() != InitialSpace->GetHilbertSpaceDimension())
    {
      cout << "error, initial Hilbert space dimension is "  << InitialSpace->GetHilbertSpaceDimension() << " while " << Manager.GetString("ground-file") << " has dimension " << InitialState.GetVectorDimension() << endl;
    }
  ComplexVector FinalState = InitialSpace->ParticleHoleSymmetrize(InitialState, *TargetSpace); 
  FinalState.WriteVector(OutputFileName);

  delete TargetSpace;
  delete InitialSpace;
  
  return 0;
}
