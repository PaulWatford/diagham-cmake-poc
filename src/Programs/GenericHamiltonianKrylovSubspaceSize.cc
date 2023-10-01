#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"
#include "Matrix/LongRationalMatrix.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Hamiltonian/FileBasedHamiltonian.h"
#include "Hamiltonian/FileBasedHermitianHamiltonian.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"

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
  OptionManager Manager ("GenericHamiltonianKrylovSubspaceSize" , "0.01");
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
  (*SystemGroup) += new  SingleIntegerOption ('\n', "skip-lines", "skip the first n-tf lines of the input file", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "data-column", "index of the column that contains the matrix elements (or their real part)", 2);
  (*SystemGroup) += new BooleanOption  ('s', "nosort-rowindices", "do not sort the row indices, assuming they are already sorted from the smallest to the largest"); 
  (*SystemGroup) += new SingleIntegerOption ('\n', "rank-step", "if non-zero, test the rank of the Krylov set at every rank-step", 0); 
  //  (*SystemGroup) += new BooleanOption  ('\n', "all-basis", "scan the Krylov subpsace size for each basis state"); 
  (*SystemGroup) += new  SingleIntegerOption ('\n', "basis-minrange", "scan the Krylov subpsace size for each basis state with at least this min index", 0); 
  (*SystemGroup) += new  SingleIntegerOption ('\n', "basis-maxrange", "scan the Krylov subpsace size for each basis state with up to this max index (if negative, up to Hilbert space dimension)", -1); 
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "output name for the characteristic polynomial");
  (*MiscGroup) += new BooleanOption  ('\n', "show-time", "show the amount of time spent in various steps of the calculation");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericHamiltonianKrylovSubspaceSize -h" << endl;
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

  int PartialRankTest = Manager.GetInteger("rank-step");

  if (false)
    {
      AbstractHamiltonian* Hamiltonian = 0;
      Hamiltonian  = new FileBasedHamiltonian(Manager.GetString("hamiltonian"), Manager.GetInteger("data-column"), false, Manager.GetBoolean("base-one"),
					      Manager.GetInteger("skip-lines"), Manager.GetBoolean("nosort-rowindices"));
      
      cout << "Hilbert space dimension = " << Hamiltonian->GetHilbertSpaceDimension() << endl;
      Architecture.GetArchitecture()->SetDimension(Hamiltonian->GetHilbertSpaceDimension());	
      
      FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture.GetArchitecture(), 1, Hamiltonian->GetHilbertSpaceDimension(), false);
      Lanczos.SetHamiltonian(Hamiltonian);
      Lanczos.InitializeLanczosAlgorithm() ;
      Lanczos.RunLanczosAlgorithm(Hamiltonian->GetHilbertSpaceDimension() - 1);
      
      return 0;
      RealMatrix Krylov (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension(), true);
      RealMatrix TmpKrylov (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension(), true);
      
      Krylov.SetMatrixElement(Hamiltonian->GetHilbertSpaceDimension() - 1, 0, 1.0);

      for (int i = 1; i < Hamiltonian->GetHilbertSpaceDimension(); ++i)
	{
	  Hamiltonian->Multiply(Krylov[i - 1], Krylov[i]);
	  Krylov[i] /= Krylov[i].Norm();
	  TmpKrylov.Copy(Krylov);
	  double* TmpValues = TmpKrylov.SingularValueDecomposition();
	  for (int j = 0; j < i; ++j)
	    {
	      cout << TmpValues[j] << endl;
	    }
	  cout << endl;
	  delete[] TmpValues;
	}
      
      double* TmpValues = Krylov.SingularValueDecomposition();
      for (int i = 0; i < Hamiltonian->GetHilbertSpaceDimension(); ++i)
	{
	  cout << TmpValues[i] << endl;
	}
      
      delete Hamiltonian;
    }
  else
    {
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
      
      //      LongIntegerMatrix TmpMatrix(NbrColumns, NbrColumns, true);
      IntegerMatrix TmpMatrix(NbrColumns, NbrColumns, true);

      long* MatrixElements = HamiltonianFile.GetAsLongArray(2);
      if (MatrixElements == 0)
	{
	  HamiltonianFile.DumpErrors(cout) << endl;
	  return 0;
	}
      for (long i = 0l; i < NbrMatrixElements; ++i)
	{
	  TmpMatrix.SetMatrixElement(RowIndices[i], ColumnIndices[i], MatrixElements[i]);
	}

      //      if (Manager.GetBoolean("all-basis"))
      //	{

      int MinRange = Manager.GetInteger("basis-minrange");
      int MaxRange = Manager.GetInteger("basis-maxrange");
      if (MaxRange < 0)
	{
	  MaxRange = NbrColumns;
	}
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      timeval StartingTime;
      timeval EndingTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
      for (int j = MinRange; j < MaxRange; ++j)
	{
	  //	  LongRationalMatrix Krylov(NbrColumns, NbrColumns, true);
	  LongIntegerMatrix Krylov(NbrColumns, NbrColumns, true);
	  Krylov.SetMatrixElement(j, 0, 1);
	  int TmpRank = 0;
	  for (int i = 1; ((i < NbrColumns) && (TmpRank == 0)); ++i)
	    {
	      Krylov[i].Multiply(TmpMatrix, Krylov[i - 1]);
	      if ((PartialRankTest > 0) && ((i % PartialRankTest) == 0))
		{
		  //		  LongRationalMatrix Krylov2(NbrColumns, NbrColumns, true);
		  //		  Krylov2.Copy(Krylov);
		  LongRationalMatrix Krylov2 (Krylov);
		  gettimeofday (&(StartingTime), 0);
		  Krylov2.Resize(NbrColumns, i + 1);
		  int TmpRank2 = Krylov2.Rank();
		  if (TmpRank2 != (i + 1))
		    {
		      TmpRank = TmpRank2;
		    }
		  gettimeofday (&(EndingTime), 0);
		  Dt = (double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
		    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1.0e6);		      
		  if (Manager.GetBoolean("show-time"))
		    {
		      cout << "partial Krylov subspace dimension evaluated in " << Dt << "s" << endl;
		    }
		}
	    }
	  if (TmpRank == 0)
	    {
	      TmpRank = Krylov.Rank();
	    }
	  cout << "size Krylov subspace |" << j << "> = " << TmpRank << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1.0e6);		      
	  if (Manager.GetBoolean("show-time"))
	    {
	      cout << "Krylov subspace dimension evaluated in " << Dt << "s" << endl;
	    }
	}
      //    }
    }
  return 0;
}
