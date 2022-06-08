#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"

#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslationsAndSublatticeConservation.h"

#include "Hamiltonian/ParticleOnTorusPairHoppingHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusPairHoppingRealHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"
#include "MainTask/FQHEOnTorusMainTask.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstdio>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
void FQHPairHoppingComputeCharacteristicPolynomial(AbstractQHEHamiltonian* hamiltonian, ParticleOnTorusWithMagneticTranslations* chain, char* outputFileName, AbstractArchitecture* architecture);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHETorusFermionsPairHopping" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-sites", "number of sites", 18);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "file describing the interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the full reduced Brillouin zone", false);
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new  BooleanOption ('\n', "enable-realhamiltonian", "use a real Hamiltonian at the inversion symmetric points");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-sublatticeconservation", "for even number of sites, do not use the sublattice particle number conservation");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "export-charpolynomial", "export the hamiltonian characteristic polynomial");  
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsPairHopping -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("nbr-sites");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  char *LoadPrecalculationFile=Manager.GetString("load-precalculation");
  int LandauLevel=0;
  int NbrPseudopotentials=0;
  double *Pseudopotentials=NULL;
  double HaveCoulomb=false;
  char *InteractionName=NULL;

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* SuffixOutputName = new char [256];
  sprintf (SuffixOutputName, "n_%d_2s_%d.dat", NbrFermions, MaxMomentum);

  char* OutputName = new char [512 + strlen(SuffixOutputName)];
  sprintf (OutputName, "fermions_pairhopping_%s", SuffixOutputName);

  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  cout << "test " << MomentumModulo << endl;
  int XMaxMomentum = (MomentumModulo - 1);
  bool GenerateMomenta = false;
  if ((XMomentum < 0)||(YMomentum < 0))
    GenerateMomenta = true;
  if (XMomentum < 0)
    XMomentum = 0;
  else
    XMaxMomentum = XMomentum;
  int YMaxMomentum = (MaxMomentum - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum;

  int NbrMomenta;
  int *XMomenta;
  int *YMomenta;
  int *Multiplicities = NULL;
  int CenterX=0, CenterY=0;

  if (GenerateMomenta == false)
    {
      NbrMomenta=1;
      XMomenta = new int[1];
      YMomenta = new int[1];
      XMomenta[0] = XMomentum;
      YMomenta[0] = YMomentum;
    }
  else
    {
      if (Manager.GetBoolean("all-points"))
	{
	  int Pos=0;
	  NbrMomenta = (XMaxMomentum-XMomentum+1)*(YMaxMomentum-YMomentum+1);
	  XMomenta = new int[NbrMomenta];
	  YMomenta = new int[NbrMomenta];
	  for (; XMomentum <= XMaxMomentum; ++XMomentum)
	    for (int YMomentum2 = YMomentum; YMomentum2<= YMaxMomentum; ++YMomentum2)
	      {
		XMomenta[Pos]=XMomentum;
		YMomenta[Pos]=YMomentum2;
		++Pos;
	      }
	}
      else // determine inequivalent states in BZ
	{
	  if (Manager.GetBoolean("full-reducedbz"))
	    {
	      int Pos=0;
	      XMaxMomentum = MomentumModulo;
	      YMaxMomentum = MomentumModulo;
	      NbrMomenta = MomentumModulo * MomentumModulo;
	      XMomenta = new int[NbrMomenta];
	      YMomenta = new int[NbrMomenta];
	      for (; XMomentum < XMaxMomentum; ++XMomentum)
		{
		  for (int YMomentum2 = YMomentum; YMomentum2 < YMaxMomentum; ++YMomentum2)
		    {
		      XMomenta[Pos] = XMomentum;
		      YMomenta[Pos] = YMomentum2;
		      ++Pos;
		    }
		}
	    }
	  else
	    {
	      if (NbrFermions & 1)
		{
		  CenterX=0;
		  CenterY=0;
		}
	      else
		{
		  if ((NbrFermions/MomentumModulo*MaxMomentum/MomentumModulo)&1) // p*q odd?
		    {
		      CenterX=MomentumModulo/2;
		      CenterY=MomentumModulo/2;
		    }
		  else
		    {
		      CenterX=0;
		      CenterY=0;
		    }
		}
	      NbrMomenta=0;
	      for (int Kx = CenterX; Kx<=CenterX+MomentumModulo/2; ++Kx)
		for (int Ky= (Kx-CenterX)+CenterY; Ky<=CenterY+MomentumModulo/2; ++Ky)
		  {
		    ++NbrMomenta;
		  }
	      int Pos=0;
	      XMomenta = new int[NbrMomenta];
	      YMomenta = new int[NbrMomenta];
	      Multiplicities = new int[NbrMomenta];
	      for (int Kx = 0; Kx<=MomentumModulo/2; ++Kx)
		for (int Ky= Kx; Ky<=MomentumModulo/2; ++Ky, ++Pos)
		  {
		    XMomenta[Pos]=CenterX+Kx;
		    YMomenta[Pos]=CenterY+Ky;
		    if (Kx==0)
		      {
			if (Ky==0)
			  Multiplicities[Pos]=1; // BZ center
			else if (Ky==MomentumModulo/2)
			  Multiplicities[Pos]=2;
			else Multiplicities[Pos]=4;
		      }
		    else if (Kx==MomentumModulo/2)
		      {
			Multiplicities[Pos]=1; // BZ corner
		      }
		    else
		      {
			if (Ky==Kx) // diagonal ?
			  {
			    Multiplicities[Pos]=4; 
			  }
			else
			  {
			    if (Ky==MomentumModulo/2)
			      Multiplicities[Pos]=4;
			    else
			      Multiplicities[Pos]=8;
			  }
		      }
		  }	    
	    }
	}
    }
  
  char* CommentLine = new char [512];
  if (((MaxMomentum & 1) == 0) && (Manager.GetBoolean("disable-sublatticeconservation") == false))
    {
      sprintf (CommentLine, " periodic pair hopping chain with %d sites \n# Kx Ky Sub Energy ", MaxMomentum);
    }
  else
    {
       sprintf (CommentLine, " periodic pair hopping chain with %d sites \n# Kx Ky Energy ", MaxMomentum);
   }
  
  bool FirstRun = true;
  bool ForceRealFlag = false;
  for (int Pos = 0; Pos < NbrMomenta; ++Pos)
    {
      XMomentum = XMomenta[Pos];
      YMomentum = YMomenta[Pos];
      int MaxNbrFermionsEvenMomentum = 0;
      if (((MaxMomentum & 1) == 0) && (Manager.GetBoolean("disable-sublatticeconservation") == false))
	{
	  MaxNbrFermionsEvenMomentum = NbrFermions;
	}
      for (int NbrFermionsEvenMomentum = 0; NbrFermionsEvenMomentum <= MaxNbrFermionsEvenMomentum; ++NbrFermionsEvenMomentum)
	{
	  cout << "----------------------------------------------------------------" << endl;
	  cout << "kx=" << XMomentum << ", ky=" << YMomentum;
	  if (MaxNbrFermionsEvenMomentum != 0)
	    {
	      cout << ", even sub=" << NbrFermionsEvenMomentum;
	    }
	  cout << endl;
	  FermionOnTorusWithMagneticTranslations* TotalSpace = 0;
	  if (MaxNbrFermionsEvenMomentum != 0)
	    {
	      TotalSpace = new FermionOnTorusWithMagneticTranslationsAndSublatticeConservation(NbrFermions, MaxMomentum, XMomentum, YMomentum, NbrFermionsEvenMomentum);
	    }
	  else
	    {
	      TotalSpace = new FermionOnTorusWithMagneticTranslations(NbrFermions, MaxMomentum, XMomentum, YMomentum);
	    }
	  //cout << " Total Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
	  //cout << "momentum = (" << XMomentum << "," << YMomentum << ")" << endl;
	  
	  if (TotalSpace->GetHilbertSpaceDimension() > 0)
	    {
	      
	      Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());
	      
	      AbstractQHEHamiltonian* Hamiltonian;
	      if ((Manager.GetBoolean("enable-realhamiltonian") == true) &&
		  (((XMomentum % MomentumModulo) == 0) || (((MomentumModulo & 1) == 0) && ((XMomentum % MomentumModulo) == (MomentumModulo / 2)))) &&
		  (((YMomentum % MomentumModulo) == 0) || (((MomentumModulo & 1) == 0) && ((YMomentum % MomentumModulo) == (MomentumModulo / 2)))))
		{
		  cout << "using real hamiltonian" << endl;
		  ForceRealFlag = true;
		  Lanczos.SetRealAlgorithms();
		  Hamiltonian = new ParticleOnTorusPairHoppingRealHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
									      Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
		}
	      else
		{		      
		  cout << "using complex hamiltonian" << endl;
		  Hamiltonian = new ParticleOnTorusPairHoppingHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
									  Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
		}
	      
	      char* EigenvectorName = 0;
	      if ((Manager.GetBoolean("eigenstate")) || (Manager.GetBoolean("export-charpolynomial")))
		{
		  EigenvectorName = new char [512];
		  char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
		  if (MaxNbrFermionsEvenMomentum != 0)
		    {
		      if (Manager.GetString("eigenstate-file") == 0)
			sprintf (EigenvectorName, "%s_sublat_%d_kx_%d_ky_%d", TmpName, NbrFermionsEvenMomentum, XMomentum, YMomentum);
		      else
			sprintf (EigenvectorName, "%s_sublat_%d_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), NbrFermionsEvenMomentum, XMomentum, YMomentum);
		    }
		  else
		    {
		      if (Manager.GetString("eigenstate-file") == 0)
			sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
		      else
			sprintf (EigenvectorName, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), XMomentum, YMomentum);
		    }
		  delete [] TmpName;
		}
	      if (Manager.GetBoolean("export-charpolynomial"))
		{
		  FQHPairHoppingComputeCharacteristicPolynomial(Hamiltonian, TotalSpace, EigenvectorName, Architecture.GetArchitecture());
		}
	      
	      double Shift = 0.0;
	      Hamiltonian->ShiftHamiltonian(Shift);
	      char* TmpSzString = new char[64];
	      if (MaxNbrFermionsEvenMomentum != 0)
		{
		  sprintf (TmpSzString, "%d %d %d", XMomentum, YMomentum, NbrFermionsEvenMomentum);
		}
	      else
		{
		  sprintf (TmpSzString, "%d %d", XMomentum, YMomentum);
		}
	      if (ForceRealFlag == true)
		{
		  GenericRealMainTask Task(&Manager, TotalSpace, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  OutputName,
					   FirstRun, EigenvectorName);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		}
	      else
		{
		  GenericComplexMainTask Task(&Manager, TotalSpace, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  OutputName,
					      FirstRun, EigenvectorName);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		}
	      delete[] TmpSzString;
	      //	      FQHEOnTorusMainTask Task (&Manager, TotalSpace, &Lanczos, Hamiltonian, YMomentum, Shift, OutputName, FirstRun, EigenvectorName, XMomentum, 0, ForceRealFlag);
	      //      Task.SetKxValue(XMomentum);
	      // if (Multiplicities!=0)
	      // 	Task.SetMultiplicity(Multiplicities[Pos]);
	      if (EigenvectorName != 0)
		{
		  delete[] EigenvectorName;
		}
	      if (FirstRun == true)
		FirstRun = false;
	      Lanczos.SetComplexAlgorithms();
	      ForceRealFlag = false;
	      delete Hamiltonian;
	    }
	  delete TotalSpace;
	}
    }
  delete [] XMomenta;
  delete [] YMomenta;
  if (Multiplicities!=0)
    delete [] Multiplicities;
  return 0;
}


// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// discardFourFactor = discard a global four factor used to ensure integer numbers

void FQHPairHoppingComputeCharacteristicPolynomial(AbstractQHEHamiltonian* hamiltonian, ParticleOnTorusWithMagneticTranslations* chain, char* outputFileName, AbstractArchitecture* architecture)
{
#ifdef __GMP__
  cout << "Computing the hamiltonian" << endl;
  RealMatrix TmpRawMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
  hamiltonian->GetHamiltonian(TmpRawMatrix);
  double* TmpNormalizationFactors = chain->GetBasisNormalization();
  for (int i = 0; i < chain->GetHilbertSpaceDimension(); ++i)
    {
      for (int j = 0; j < chain->GetHilbertSpaceDimension(); ++j)
	{
	  double Tmp;
	  TmpRawMatrix.GetMatrixElement(i, j, Tmp);
	  Tmp *= TmpNormalizationFactors[i];
	  Tmp /= TmpNormalizationFactors[j];
	  TmpRawMatrix.SetMatrixElement(i, j, Tmp);
	}
    }
  cout << "Converting to integer matrix" << endl;
  LongIntegerMatrix TmpMatrix;
  TmpMatrix = LongIntegerMatrix(TmpRawMatrix);

  cout << "Start computing characteristic polynomial (degree " << chain->GetHilbertSpaceDimension() << ")" << endl;
  mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(architecture);
  char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
  sprintf (PolynomialOutputFileName, "%s.charpol", outputFileName);
  ofstream OutputFile;
  OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
  OutputFile << CharacteristicPolynomial[0];
  for (int i = 1; i <= chain->GetHilbertSpaceDimension(); ++i)
    {
      OutputFile << "," << CharacteristicPolynomial[i];
    }
  OutputFile << endl;
  OutputFile.close();
#else
  cout << "GMP library is required for characteristic polynomials" << endl;
#endif	       
}

