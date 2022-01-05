#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"
#include "Hamiltonian/SpinChainJz2Hamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"

#include "HilbertSpace/Spin3_2Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainWithIntegerHamiltonian" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "use-mirror", "use the mirror symmetry");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainWithIntegerHamiltonian -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "integer_spin_%d_%schain_n_%d", (SpinValue / 2), BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-mirror") == true)
	{
	  sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz P_mirror ", BoundaryName, (SpinValue / 2), NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz ", BoundaryName, (SpinValue / 2), NbrSpins);
	}
    }
  else
    {
      sprintf (OutputFileName, "integer_spin_%d_2_%schain_n_%d", SpinValue, BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-mirror") == true)
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz P_mirror", BoundaryName, SpinValue, NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz", BoundaryName, SpinValue, NbrSpins);
	}
    }
  char* FullOutputFileName = new char [strlen(OutputFileName) + 64];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  double* JValues = new double [NbrSpins];
  if ((SpinValue & 1) == 0)
    {
      JValues[0] = 1.0;
    }
  else
    {
      JValues[0] = 4.0;
    }
  for (int i = 1; i < NbrSpins; ++i)
    JValues[i] = JValues[0];
  double* JzValues = new double [NbrSpins];
  for (int i = 0; i < NbrSpins; ++i)
    JzValues[i] = JValues[i];
  double* HzValues = 0;

  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if (Manager.GetInteger("initial-sz") != 0)
    {
      InitalSzValue = Manager.GetInteger("initial-sz");
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;
  if (Manager.GetBoolean("use-mirror") == true)
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  for (int Mirror = 0; Mirror < 2; ++Mirror)
	    { 
	      AbstractSpinChain* Chain = 0;
	      switch (SpinValue)
		{
		case 1:
		  {
		    Chain = new Spin1_2ChainMirrorSymmetry (NbrSpins, InitalSzValue, Mirror, 1000000);
		  }
		  break;
		case 2:
		  {
		    Chain = new Spin1ChainWithInversionSymmetry (NbrSpins, (2 * Mirror) - 1, InitalSzValue, 1000000);
		  }
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{
		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		  SpinChainHamiltonian* Hamiltonian = 0;
		  if (HzValues == 0)
		    Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, Manager.GetBoolean("use-periodic"));
		  else
		    Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, HzValues, Manager.GetBoolean("use-periodic"));
		  char* TmpSzString = new char[64];
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpSzString, "%d %d", InitalSzValue, Mirror);
		  sprintf (TmpEigenstateString, "%s_sz_%d_m_%d", OutputFileName, InitalSzValue, Mirror);


		  // IntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		  // Hamiltonian->GetHamiltonian(TmpMatrix);
		  // long* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();

		  LongIntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		  Hamiltonian->GetHamiltonian(TmpMatrix);
#ifdef __GMP__
		  mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();
		  char* PolynomialOutputFileName = new char[strlen(OutputFileName) + 256];
		  sprintf (PolynomialOutputFileName, "%s_sz_%d_m_%d.charpol", OutputFileName, InitalSzValue, Mirror);		  
		  ofstream OutputFile;
		  OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
		  OutputFile << CharacteristicPolynomial[0];
		  for (int i = 1; i <= Chain->GetHilbertSpaceDimension(); ++i)
		    {
		      OutputFile << "," << CharacteristicPolynomial[i];
		    }
		  OutputFile << endl;
		  OutputFile.close();
#else
		  cout << "GMP library is required" << endl;
#endif		 


		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete Hamiltonian;
		  delete[] TmpSzString;
		  delete[] TmpEigenstateString;
		}
	      delete Chain;
	    }
	}
    }
  else
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  AbstractSpinChain* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      {
		Chain = new Spin1_2Chain (NbrSpins, InitalSzValue, 1000000);
	      }
	      break;
	    case 2 :
	      Chain = new Spin1Chain (NbrSpins, InitalSzValue, 1000000);
	      break;
	    case 3 :
	      Chain = new Spin3_2Chain (NbrSpins, InitalSzValue, 1000000);
	      break;
	    default :
	      {
		if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
	    else 
	      cout << "spin " << SpinValue << "/2 are not available" << endl;
		return -1;
	      }
	    }
	  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());		  
	  SpinChainHamiltonian* Hamiltonian = 0;
	  if (HzValues == 0)
	    {
	      Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, Manager.GetBoolean("use-periodic"));
	    }
	  else
	    {
	      Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, HzValues, Manager.GetBoolean("use-periodic"));
	    }
	  char* TmpSzString = new char[64];
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpSzString, "%d", InitalSzValue);
	  sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);

	  
	  // IntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	  // Hamiltonian->GetHamiltonian(TmpMatrix);
	  // long* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();

	  
	  LongIntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	  Hamiltonian->GetHamiltonian(TmpMatrix);
#ifdef __GMP__
	  mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();


	  char* PolynomialOutputFileName = new char[strlen(OutputFileName) + 256];
	  sprintf (PolynomialOutputFileName, "%s_sz_%d.charpol", OutputFileName, InitalSzValue);
	  ofstream OutputFile;
	  OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
	  OutputFile << CharacteristicPolynomial[0];
	  for (int i = 1; i <= Chain->GetHilbertSpaceDimension(); ++i)
	    {
	      OutputFile << "," << CharacteristicPolynomial[i];
	    }
	  OutputFile << endl;
	  OutputFile.close();
#else
		  cout << "GMP library is required" << endl;
#endif		 
	  	    
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName, FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Hamiltonian;
	  delete Chain;
	  delete[] TmpSzString;
	  delete[] TmpEigenstateString;
	}
    }
  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] JValues;
  delete[] FullOutputFileName;
  delete[] JzValues;
  return 0;
}
