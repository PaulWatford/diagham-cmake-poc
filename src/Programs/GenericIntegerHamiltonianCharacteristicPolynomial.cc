#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/FilenameTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("GenericIntegerHamiltonianCharacteristicPolynomial" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('\n', "hamiltonian", "text file where the hamiltonian matrix elements are stored");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "shift", "shift the hamiltonian by a constant integer value (i.e., H-lambda 1)", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "base-one", "hamiltonian indices start at one rather than zero");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "output name for the characteristic polynomial");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericIntegerHamiltonianCharacteristicPolynomial -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("hamiltonian") == 0)
    {
      cout << "no hamiltonian provided" << endl; 
      return -1;
    }
  if (IsFile(Manager.GetString("hamiltonian")) == false)
    {
      cout << "can't open hamiltonian " << Manager.GetString("hamiltonian") << endl; 
      return -1;
    }

  MultiColumnASCIIFile HamiltonianFile;
  if (HamiltonianFile.Parse(Manager.GetString("hamiltonian")) == false)
    {
      return -1;
    }
  if (HamiltonianFile.GetNbrColumns() != 3)
    {
      cout << Manager.GetString("hamiltonian") << " has a wrong number of columns" << endl;
      return -1;
    }

  long NbrMatrixElements = HamiltonianFile.GetNbrLines();
    
  int* RowIndices = HamiltonianFile.GetAsIntegerArray(0);
  int* ColumnIndices = HamiltonianFile.GetAsIntegerArray(1);
  
  int NbrColumns = 0;
  if (Manager.GetBoolean("base-one"))
    {
      for (long i = 0l; i < NbrMatrixElements; ++i)
	{
	  if (RowIndices[i] > NbrColumns)
	    {
	      NbrColumns = RowIndices[i];
	    }
	  if (RowIndices[i] == 0)
	    {
	      cout << "error on entry " << i << " (" << RowIndices[i] << " " << ColumnIndices[i] << "), zero index found while using --base-one option" << endl;
	      return 0;
	    }
	  RowIndices[i]--;
	  if (ColumnIndices[i] > NbrColumns)
	    {
	      NbrColumns = ColumnIndices[i];
	    }
	  if (ColumnIndices[i] == 0)
	    {
	      cout << "error on entry " << i << " (" << RowIndices[i] << " " << ColumnIndices[i] << "), zero index found while using --base-one option" << endl;
	      return 0;
	    }
	  ColumnIndices[i]--;
	}
    }
  else
    {
      for (long i = 0l; i < NbrMatrixElements; ++i)
	{
	  if (RowIndices[i] > NbrColumns)
	    {
	      NbrColumns = RowIndices[i];
	    }
	  if (ColumnIndices[i] > NbrColumns)
	    {
	      NbrColumns = ColumnIndices[i];
	    }
	}
      NbrColumns++;
    }
  cout << "Matrix size = " << NbrColumns << "x" << NbrColumns << " (" << NbrMatrixElements << " nbr matrix elements)" << endl;
  
  LongIntegerMatrix TmpMatrix(NbrColumns, NbrColumns, true);

  long* MatrixElements = HamiltonianFile.GetAsLongArray(2);
  for (long i = 0l; i < NbrMatrixElements; ++i)
    {
      TmpMatrix.SetMatrixElement(RowIndices[i], ColumnIndices[i], MatrixElements[i]);
    }

  long EnergyShift = Manager.GetInteger("shift");
  if (EnergyShift != 0l)
    {
      for (long i = 0l; i < NbrColumns; ++i)
	{
	  TmpMatrix.AddToMatrixElement(i, i, -EnergyShift);
	}
    }

  char* PolynomialOutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      if (EnergyShift != 0l)
	{
	  PolynomialOutputFileName = new char[strlen(Manager.GetString("hamiltonian")) + 64];
	  sprintf (PolynomialOutputFileName, "%s_shift_%ld.charpol", Manager.GetString("hamiltonian"), EnergyShift);
	}
      else
	{
	  PolynomialOutputFileName = new char[strlen(Manager.GetString("hamiltonian")) + 16];
	  sprintf (PolynomialOutputFileName, "%s.charpol", Manager.GetString("hamiltonian"));
	}
    }
  else
    {
      PolynomialOutputFileName = new char[strlen(Manager.GetString("output-file")) + 16];
      strcpy (PolynomialOutputFileName, Manager.GetString("output-file"));
    }

  Architecture.GetArchitecture()->SetDimension(NbrColumns);

  if (Architecture.GetArchitecture()->CanWriteOnDisk())
    {
#ifdef __GMP__      
      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(Architecture.GetArchitecture());
      ofstream OutputFile;
      OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
      
      OutputFile << CharacteristicPolynomial[0];
      for (int i = 1; i <= NbrColumns; ++i)
	{
	  OutputFile << "," << CharacteristicPolynomial[i];
	}
      OutputFile << endl;
      OutputFile.close();
#else
      cout << "GMP library is required" << endl;
#endif      
    }
  delete[] PolynomialOutputFileName;
  
  return 0;
}
