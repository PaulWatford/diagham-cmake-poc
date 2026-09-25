#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"
#include "Matrix/LongRationalMatrix.h"

#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslationsAndSublatticeConservation.h"

#include "Hamiltonian/ParticleOnTorusPairHoppingHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusPairHoppingRealHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
#include <limits>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


// compute the characteristic polynomial for hamiltonian projected onto the Krylov subspace
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// productStateIndices = indices of the product states that defines the Krylov subpspace
// productStateCoefficients = integer coefficients for each product state of the Krylov subpspace generating vector
// nbrProductStateConfigurations = number o product states defining the Krylov subpspace generating vector
// exportKrylov = if true, export the krylov subspace in a tex file
// exportMathematicaKrylov = if true, export the krylov subspace in a mathematica friendly text file
// return value = Krylov subspace dimension
int FQHPairHoppingKrylovSubspaceComputeCharacteristic(AbstractQHEHamiltonian* hamiltonian, ParticleOnTorusWithMagneticTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, int* productStateIndices, long* productStateCoefficients, int nbrProductStateConfigurations, bool exportKrylov, bool exportMathematicaKrylov);

// convert a string describing a product state in the (u,d,+,-) basis to its fermionic occupation version
//
// productState = product state in the (u,d,+,-) basis
// nbrFermions= reference to the number of fermions in the product state
// nbrFermionsEvenMomentum = reference to the sublattice sector
// maxMomentum = reference to the max momentumw number or number of sites
// yMomentum = reference to the y momentumw
// return value = array containing the fermionic occupation (0 if an error occurred) 
int* FQHPairHoppingKrylovParseProductState(char* productState, int& nbrFermions, int& nbrFermionsEvenMomentum, int& maxMomentum, int& yMomentum);


int main(int argc, char** argv)
{
  cout.precision(std::numeric_limits<double>::max_digits10);

  // some running options and help
  OptionManager Manager ("FQHETorusFermionsPairHoppingKrylovSubspace" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('\n', "product-state", "string describing the product state defining the Krylov subspace");
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-productstates", "file describing the state generating the Krylov subspace");  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-sites", "number of sites", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction", 0);
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new  BooleanOption ('\n', "enable-realhamiltonian", "use a real Hamiltonian at the inversion symmetric points");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-sublatticeconservation", "for even number of sites, do not use the sublattice particle number conservation");
  (*SystemGroup) += new  BooleanOption ('\n', "use-double", "also perform the Krylov projected calculation using double precision");
  (*SystemGroup) += new  BooleanOption ('\n', "compute-eigenstates", "when using double precision, compute the eigenstates within the Krylov subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "export-krylov", "export the krylov subspace in a tex file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-mathematica", "export the krylov subspace in a mathemtica text file");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsPairHoppingKrylovSubspace -h" << endl;
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
  int NbrFermionsEvenMomentum = 0;
  char* LoadPrecalculationFile = Manager.GetString("load-precalculation");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if ((Manager.GetString("product-state") == 0) && (Manager.GetString("multiple-productstates") == 0))
    {
      cout << "error, a product state should be provided" << endl;
      return 0;
    }
  if ((Manager.GetString("product-state") != 0) && (Manager.GetString("multiple-productstates") != 0))
    {
      cout << "error, options --product-state and --multiple-productstates cannot be used simultaneously" << endl;
      return 0;
    }

  int** ProductStateConfigurations = 0;
  int NbrProductStateConfigurations = 0;
  long* ProductStateCoefficients = 0;
  char** ProductStateConfigurationNames = 0;
  if (Manager.GetString("product-state") != 0)
    {
      NbrProductStateConfigurations = 1;
      ProductStateConfigurations = new int* [NbrProductStateConfigurations];
      ProductStateConfigurations[0] = FQHPairHoppingKrylovParseProductState(Manager.GetString("product-state"), NbrFermions, NbrFermionsEvenMomentum, MaxMomentum, YMomentum); 
      if (ProductStateConfigurations[0] == 0)
	{
	  return 0;
	}
      ProductStateConfigurationNames = new char*[NbrProductStateConfigurations];
      ProductStateConfigurationNames[0] = new char[strlen(Manager.GetString("product-state")) + 1];
      strcpy (ProductStateConfigurationNames[0], Manager.GetString("product-state"));
      ProductStateCoefficients = new long [NbrProductStateConfigurations];
      ProductStateCoefficients[0] = 1l;
    }
  else
    {
      MultiColumnASCIIFile InputFile(':');
      if (InputFile.Parse(Manager.GetString("multiple-productstates")) == false)
	{
	  InputFile.DumpErrors(cout) << endl;
	  return -1;
	}
      NbrProductStateConfigurations = InputFile.GetNbrLines();
      ProductStateConfigurationNames = new char*[NbrProductStateConfigurations];
      ProductStateCoefficients = InputFile.GetAsLongArray(0);
      ProductStateConfigurations = new int* [NbrProductStateConfigurations];
      
      for (int i = 0; i < NbrProductStateConfigurations; ++i)
	{
	  int TmpNbrFermions = 0;
	  int TmpNbrFermionsEvenMomentum = 0;
	  int TmpMaxMomentum = 0;
	  int TmpYMomentum = 0;
	  ProductStateConfigurations[i] = FQHPairHoppingKrylovParseProductState(InputFile(1, i), TmpNbrFermions,
										TmpNbrFermionsEvenMomentum, TmpMaxMomentum, TmpYMomentum); 
	  if (ProductStateConfigurations[i] == 0)
	    {
	      return 0;
	    }
	  if (i == 0)
	    {
	      NbrFermions = TmpNbrFermions;
	      NbrFermionsEvenMomentum = TmpNbrFermionsEvenMomentum;
	      MaxMomentum = TmpMaxMomentum;
	      YMomentum = TmpYMomentum;
	    }
	  else
	    {
	      if (NbrFermions != TmpNbrFermions)
		{
		  cout << InputFile(1, 0) << " and " << InputFile(1, i) << " have a different number of fermions" << endl;
		  return 0;
		}
	      if (MaxMomentum != TmpMaxMomentum)
		{
		  cout << InputFile(1, 0) << " and " << InputFile(1, i) << " have a different number of sites" << endl;
		  return 0;
		}
	      if (YMomentum != TmpYMomentum)
		{
		  cout << InputFile(1, 0) << " and " << InputFile(1, i) << " have a different momentum" << endl;
		  return 0;
		}
	      if (NbrFermionsEvenMomentum != TmpNbrFermionsEvenMomentum)
		{
		  cout << InputFile(1, 0) << " and " << InputFile(1, i) << " are in a different sublattice sector" << endl;
		  return 0;
		}
	    }
	  ProductStateConfigurationNames[i] = new char[strlen(InputFile(1, i)) + 1];
	  strcpy (ProductStateConfigurationNames[i], InputFile(1, i));	  
	}
    }
  
  char* SuffixOutputName = new char [256];
  sprintf (SuffixOutputName, "n_%d_2s_%d.dat", NbrFermions, MaxMomentum);

  char* OutputName = new char [512 + strlen(SuffixOutputName)];
  if (NbrProductStateConfigurations == 1)
    {
      sprintf (OutputName, "fermions_pairhopping_krylov_%s_%s", ProductStateConfigurationNames[0], SuffixOutputName);
    }
  else
    {
      sprintf (OutputName, "fermions_pairhopping_krylov_multiple_%s_%s", ProductStateConfigurationNames[0], SuffixOutputName);
    }
  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  int XMaxMomentum = (MomentumModulo - 1);
  bool GenerateMomenta = false;
  if (XMomentum < 0)
    XMomentum = 0;
  if (YMomentum < 0)
    YMomentum = 0;
  if ((XMomentum != 0) && (((MomentumModulo & 1) != 0) || ((XMomentum % MomentumModulo) != (MomentumModulo / 2))))
    {
      cout << "error, only kx=0 or pi are valid" << endl;
    }
  if ((YMomentum != 0) && (((MomentumModulo & 1) != 0) || ((YMomentum % MomentumModulo) != (MomentumModulo / 2))))
    {
      cout << "error, only ky=0 or pi are valid" << endl;
      return 0;
    }
  
  int NbrMomenta = 1;
  int* XMomenta = new int[1];
  int* YMomenta = new int[1];

  XMomenta[0] = XMomentum;
  YMomenta[0] = YMomentum;
  
  char* CommentLine = new char [1024];
  if (((MaxMomentum & 1) == 0) && (Manager.GetBoolean("disable-sublatticeconservation") == false))
    {
      sprintf (CommentLine, " periodic pair hopping chain with %d sites\n# for the Krylov subspace generated by \n# Kx Ky Sub Energy ", MaxMomentum);
    }
  else
    {
       sprintf (CommentLine, " periodic pair hopping chain with %d sites\n# for the Krylov subspace generated by \n# Kx Ky Energy ", MaxMomentum);
   }
  
  bool FirstRun = true;
  for (int Pos = 0; Pos < NbrMomenta; ++Pos)
    {
      XMomentum=XMomenta[Pos];
      YMomentum=YMomenta[Pos];
      int MaxNbrFermionsEvenMomentum = 0;
      if (((MaxMomentum & 1) == 0) && (Manager.GetBoolean("disable-sublatticeconservation") == false))
	{
	  MaxNbrFermionsEvenMomentum = NbrFermions;
	}
      FermionOnTorusWithMagneticTranslations* TotalSpace = 0;
      if (MaxNbrFermionsEvenMomentum != 0)
	{
	  TotalSpace = new FermionOnTorusWithMagneticTranslationsAndSublatticeConservation(NbrFermions, MaxMomentum, XMomentum, YMomentum, NbrFermionsEvenMomentum);
	}
      else
	{
	  TotalSpace = new FermionOnTorusWithMagneticTranslations(NbrFermions, MaxMomentum, XMomentum, YMomentum);
	}
      
      if (TotalSpace->GetHilbertSpaceDimension() > 0)
	{
	  int* ProductStateIndices = new int[NbrProductStateConfigurations];
	  for (int i = 0; i < NbrProductStateConfigurations; ++i)
	    {
	      ProductStateIndices[i] = TotalSpace->FindStateIndex(ProductStateConfigurations[i]);
	      if (ProductStateIndices[i] == -1)
		{
		  cout << "error, the product state configuration \"" << ProductStateConfigurationNames[i] << "\" is not compatible with the Hilbert space" << endl;
		  return 0;
		}
	      cout << "Product state configuration " << ProductStateConfigurationNames[i] << " converted as ";
	      TotalSpace->PrintState(cout, ProductStateIndices[i]) << " (index=" << ProductStateIndices[i] << ", coefficient=" << ProductStateCoefficients[i] << ")" << endl;
	    }
	  Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());
	  
	  AbstractQHEHamiltonian* Hamiltonian;
	  Lanczos.SetRealAlgorithms();
	  Hamiltonian = new ParticleOnTorusPairHoppingRealHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
								      Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
	  double Shift = 0.0;
	  Hamiltonian->ShiftHamiltonian(Shift);
	  
	  char* EigenvectorName = 0;

	  char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	  EigenvectorName = new char [strlen(TmpName) + 512];
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
	  
	  
	  int KrylovSpaceDimension = FQHPairHoppingKrylovSubspaceComputeCharacteristic(Hamiltonian, TotalSpace, EigenvectorName, Architecture.GetArchitecture(), ProductStateIndices, ProductStateCoefficients, NbrProductStateConfigurations, Manager.GetBoolean("export-krylov"), Manager.GetBoolean("export-mathematica"));

	  if (Manager.GetBoolean("use-double"))
	    {
	      RealVector* TmpVectors = new RealVector[TotalSpace->GetHilbertSpaceDimension()];
	      TmpVectors[0] = RealVector(TotalSpace->GetHilbertSpaceDimension(), true);
	      for (int i = 0; i < NbrProductStateConfigurations; ++i)
		{
		  TmpVectors[0][ProductStateIndices[i]] = (double) ProductStateCoefficients[i];
		}
	      TmpVectors[0] /= TmpVectors[0].Norm();
	      cout << TmpVectors[0] << endl;
	      if (KrylovSpaceDimension > 1)
		{
		  FullReorthogonalizedLanczosAlgorithm TmpLanczos (Architecture.GetArchitecture(), KrylovSpaceDimension, KrylovSpaceDimension + 1);
		  TmpLanczos.SetHamiltonian(Hamiltonian);
		  TmpLanczos.InitializeLanczosAlgorithm(TmpVectors[0]);
		  TmpLanczos.RunLanczosAlgorithm(KrylovSpaceDimension);
		  RealVector* KrylovLanczosVectors = (RealVector*) TmpLanczos.GetKrylovSubspace();
		  RealMatrix KrylovLanczosVectorMatrix(KrylovLanczosVectors, KrylovSpaceDimension, true);
		  //		  cout << KrylovLanczosVectorMatrix << endl;
		  RealTriDiagonalSymmetricMatrix TmpTridiagHamiltonian = TmpLanczos.GetTridiagonalMatrix();
		  RealMatrix TridiagHamiltonian (TmpTridiagHamiltonian);
		  double* TridiagCharacteristicPolynomial = TridiagHamiltonian.CharacteristicPolynomial();
		  for (int i = 0; i <= KrylovSpaceDimension; ++i)
		    {
		      cout << TridiagCharacteristicPolynomial[i] << "x^" << i << " + "  << endl;
		    }
		  if (Manager.GetBoolean("compute-eigenstates"))
		    {
		      RealMatrix TmpEigenvectors (KrylovSpaceDimension, KrylovSpaceDimension);
		      TmpEigenvectors.SetToIdentity();
		      TmpTridiagHamiltonian.Diagonalize(TmpEigenvectors);
		      TmpTridiagHamiltonian.SortMatrixUpOrder(TmpEigenvectors);
		      RealMatrix KrylovLanczosVectorMatrix2;
		      KrylovLanczosVectorMatrix2 = KrylovLanczosVectorMatrix * TmpEigenvectors;
		      cout << "eigenvalues:" << endl;
		      RealVector TmpVector (TotalSpace->GetHilbertSpaceDimension(), true);
		      char* EigenvectorName2 = new char [strlen(EigenvectorName) + 512];
		      for (int i = 0; i < KrylovSpaceDimension; ++i)
			{
			  Hamiltonian->LowLevelMultiply(KrylovLanczosVectorMatrix2[i], TmpVector);
			  cout << TmpTridiagHamiltonian.DiagonalElement(i) << " <psi|H|psi>="
			       <<(KrylovLanczosVectorMatrix2[i] * TmpVector) << endl;
			  sprintf (EigenvectorName2, "%s.%d.vec", EigenvectorName, i);
			  KrylovLanczosVectorMatrix2[i].WriteVector(EigenvectorName2);
			}
		      delete[] EigenvectorName2;
		    }
		}
	      else
		{
		  RealVector TmpVector (TotalSpace->GetHilbertSpaceDimension(), true);
		  char* EigenvectorName2 = new char [strlen(EigenvectorName) + 512];
		  Hamiltonian->LowLevelMultiply(TmpVectors[0], TmpVector);
		  cout << "<psi|H|psi>="<< (TmpVectors[0] * TmpVector) << endl;
		  sprintf (EigenvectorName2, "%s.%d.vec", EigenvectorName, 0);
		  TmpVectors[0].WriteVector(EigenvectorName2);
		  delete[] EigenvectorName2;
		}
	    }
	  
	  delete Hamiltonian;
	}
      delete TotalSpace;
    }
  delete [] XMomenta;
  delete [] YMomenta;
  return 0;
}

// compute the characteristic polynomial for hamiltonian projected onto the Krylov subspace
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// productStateIndices = indices of the product states that defines the Krylov subpspace
// productStateCoefficients = integer coefficients for each product state of the Krylov subpspace generating vector
// nbrProductStateConfigurations = number o product states defining the Krylov subpspace generating vector
// exportKrylov = if true, export the krylov subspace in a tex file
// exportMathematicaKrylov = if true, export the krylov subspace in a mathematica friendly text file
// return value = Krylov subspace dimension

int FQHPairHoppingKrylovSubspaceComputeCharacteristic(AbstractQHEHamiltonian* hamiltonian, ParticleOnTorusWithMagneticTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, int* productStateIndices, long* productStateCoefficients, int nbrProductStateConfigurations, bool exportKrylov, bool exportMathematicaKrylov)
{
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
  LongRationalMatrix RationalHamiltonian(TmpMatrix);

  cout << "Building Krylov subspace" << endl;
  LongIntegerVector* TmpVectors = new LongIntegerVector[chain->GetHilbertSpaceDimension()];
  TmpVectors[0] = LongIntegerVector(chain->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < nbrProductStateConfigurations; ++i)
    {
#ifdef __GMP__
      mpz_set_si(TmpVectors[0][productStateIndices[i]], productStateCoefficients[i]);
#else
      TmpVectors[0][productStateIndices[i]] = productStateCoefficients[i];  
#endif
    }
  //  cout << TmpVectors[0] << endl;
  int KrylovDimension = 0;
  for (int i = 1; (i < chain->GetHilbertSpaceDimension()) && (KrylovDimension == 0); ++i)
    {
      TmpVectors[i] = LongIntegerVector(chain->GetHilbertSpaceDimension(), true);
      TmpVectors[i].Multiply(TmpMatrix, TmpVectors[i - 1]);
      if (TmpVectors[i].IsNullVector() == true)
	{
	  KrylovDimension = i;
	}
      else
	{
	  LongRationalMatrix TmpRankMatrix(TmpVectors, i + 1);
	  int TmpRank = TmpRankMatrix.Rank();
	  if (TmpRank <= i)
	    {
	      KrylovDimension = TmpRank;
	      if (exportKrylov == true)
		{
		  char* KrylovOutputFileName = new char[strlen(outputFileName) + 256];
		  sprintf (KrylovOutputFileName, "%s.krylov.tex", outputFileName);	
		  ofstream OutputFile;
		  OutputFile.open(KrylovOutputFileName, ios::binary | ios::out);
		  LongRationalMatrix TmpRankMatrix2(TmpVectors, i + 1);
		  LongRational Tmp;
		  unsigned long* TmpStateSpinBasis = new unsigned long[chain->GetNbrOrbitals()];
		  char* TmpStateSpinBasisString = new char[(chain->GetNbrOrbitals() / 2) + 2];
		  if (TmpRankMatrix2.GetNbrColumn() < 10)
		    {
		      OutputFile << "\\begin{table}[htbp] " << endl << "\\centering" << endl << "\\begin{tabular}{c|c";
		    }
		  else
		    {
		      OutputFile << "\\begin{sidewaystable}[htbp] " << endl << "\\centering" << endl << "\\begin{tabular}{c|c";
		    }
		  for (int k = 0; k < TmpRankMatrix2.GetNbrColumn(); ++k)
		    {
		      OutputFile << "|c";
		    }
		  OutputFile << "}" << endl;
		  OutputFile << "\\hline" << endl << "\\hline" << endl; 
		  for (int j = 0; j < TmpRankMatrix2.GetNbrRow(); ++j)
		    {
		      bool TmpFlag = true;
		      for (int k = 0; k < TmpRankMatrix2.GetNbrColumn(); ++k)
			{
			  if (TmpRankMatrix2[k][j].IsZero() == false)
			    {
			      TmpFlag = false;
			    }
			}
		      if (TmpFlag == false)
			{
			  for (int k = 0; k < TmpRankMatrix2.GetNbrColumn(); ++k)
			    {
			      TmpRankMatrix2.GetMatrixElement(j, k, Tmp);
			      OutputFile << "\\tiny{$" << Tmp << "$} & "; 
			    }
			  OutputFile << " \\tiny{$\\ket{";
			  chain->PrintState(OutputFile, j) << "}$}";
			  chain->GetOccupationNumber(j, TmpStateSpinBasis);
			  for (int k = 0; k < chain->GetNbrOrbitals() ; k += 2)
			    {
			      unsigned long Tmp = TmpStateSpinBasis[k] | (TmpStateSpinBasis[k + 1] << 1);			      
			      switch (Tmp)
				{
				case 0x0ul:
				  TmpStateSpinBasisString[k >> 1] = '-';
				  break;
				case 0x1ul:
				  TmpStateSpinBasisString[k >> 1] = 'd';
				  break;
				case 0x2ul:
				  TmpStateSpinBasisString[k >> 1] = 'u';
				  break;
				case 0x3ul:
				  TmpStateSpinBasisString[k >> 1] = '+';
				  break;
				}
			      TmpStateSpinBasisString[chain->GetNbrOrbitals() / 2] = '\0';
			    }
			  OutputFile << " & \\tiny{$\\ket{" << TmpStateSpinBasisString << "}$}";			  
			  OutputFile << "\\\\" << endl;
			}
		    }
		  OutputFile << "\\hline" << endl << "\\hline" << endl; 
		  OutputFile << "\\end{tabular}" << endl;
		  OutputFile << "\\caption{Krylov subspace for root=";
		  if (nbrProductStateConfigurations == 1)
		    {
		      OutputFile << "$\\ket{";
		      chain->PrintState(OutputFile, productStateIndices[0]) << "}$";
		    }
		  else
		    {
		      OutputFile << "$" << productStateCoefficients[0] << "\\ket{";
		      chain->PrintState(OutputFile, productStateIndices[0]) << "}";
		      for (int i = 0; i < nbrProductStateConfigurations; ++i)
			{
			  if (productStateCoefficients[i] < 0)
			    {
			      OutputFile << productStateCoefficients[i] << "\\ket{";
			      chain->PrintState(OutputFile, productStateIndices[0]) << "}";
			    }
			  else
			    {
			      OutputFile << "+" << productStateCoefficients[i] << "\\ket{";
			      chain->PrintState(OutputFile, productStateIndices[0]) << "}";
			    }
			}
		      OutputFile << "$";
		    }
		  OutputFile << "}" << endl;
		  if (TmpRankMatrix2.GetNbrColumn() < 10)
		    {
		      OutputFile << "\\end{table}" << endl;
		    }
		  else
		    {
		      OutputFile << "\\end{sidewaystable}" << endl;
		    }
		  OutputFile.close();
		}
	      if (exportMathematicaKrylov == true)
		{
		  char* KrylovOutputFileName = new char[strlen(outputFileName) + 256];
		  sprintf (KrylovOutputFileName, "%s.krylov.mathematica.txt", outputFileName);	
		  ofstream OutputFile;
		  OutputFile.open(KrylovOutputFileName, ios::binary | ios::out);
		  LongRationalMatrix TmpRankMatrix2(TmpVectors, i + 1);
		  LongRational Tmp;
		  unsigned long* TmpStateSpinBasis = new unsigned long[chain->GetNbrOrbitals()];
		  char* TmpStateSpinBasisString = new char[(chain->GetNbrOrbitals() / 2) + 2];
		  OutputFile << "matA={";
		  char TmpSeparator = '\0';
		  for (int j = 0; j < TmpRankMatrix2.GetNbrRow(); ++j)
		    {
		      bool TmpFlag = true;
		      for (int k = 0; k < TmpRankMatrix2.GetNbrColumn(); ++k)
			{
			  if (TmpRankMatrix2[k][j].IsZero() == false)
			    {
			      TmpFlag = false;
			    }
			}
		      if (TmpFlag == false)
			{
			  if (TmpSeparator != '\0')
			    {
			      OutputFile << TmpSeparator;
			    }
			  OutputFile << "{";
			  TmpRankMatrix2.GetMatrixElement(j, 0, Tmp);
			  OutputFile << Tmp; 
			  for (int k = 1; k < (TmpRankMatrix2.GetNbrColumn() - 1); ++k)
			    {
			      TmpRankMatrix2.GetMatrixElement(j, k, Tmp);
			      OutputFile << "," << Tmp; 
			    }
			  OutputFile << "}";
			  TmpSeparator = ',';
			}
		      else
			{
			  TmpSeparator = '\0';
			}
		    }
		  OutputFile << "}" << endl;

		  OutputFile << "vecHn={";
		  TmpSeparator = '\0';
		  for (int j = 0; j < TmpRankMatrix2.GetNbrRow(); ++j)
		    {
		      bool TmpFlag = true;
		      for (int k = 0; k < TmpRankMatrix2.GetNbrColumn(); ++k)
			{
			  if (TmpRankMatrix2[k][j].IsZero() == false)
			    {
			      TmpFlag = false;
			    }
			}
		      if (TmpFlag == false)
			{
			  if (TmpSeparator != '\0')
			    {
			      OutputFile << TmpSeparator;
			    }
			  TmpRankMatrix2.GetMatrixElement(j, TmpRankMatrix2.GetNbrColumn() - 1, Tmp);
			  OutputFile << Tmp; 
			  TmpSeparator = ',';
			}
		      else
			{
			  TmpSeparator = '\0';
			}
		    }
		  OutputFile << "}" << endl;

		  OutputFile.close();
		}
	    }
	}
    }
  cout << "Krylov subspace dimension=" << KrylovDimension << endl;

  cout << "Orthogonalizing Krylov subspace" << endl;
  LongRationalMatrix KrylovOrthogonalizeMatrix(chain->GetHilbertSpaceDimension(), KrylovDimension);
  LongRational TmpFactor;  
  KrylovOrthogonalizeMatrix[0] = TmpVectors[0];
  LongRational* NormalizationFactors = new LongRational[KrylovDimension];
  NormalizationFactors[0] = KrylovOrthogonalizeMatrix[0] * KrylovOrthogonalizeMatrix[0];
  for (int i = 1; i < KrylovDimension; ++i)
    {
        KrylovOrthogonalizeMatrix[i] = TmpVectors[i];
	for (int k = 0; k < i; ++k)
	  {
	    TmpFactor = KrylovOrthogonalizeMatrix[i] * KrylovOrthogonalizeMatrix[k];
	    TmpFactor /= NormalizationFactors[k];
	    TmpFactor.Neg();
	    KrylovOrthogonalizeMatrix[i].AddLinearCombination(TmpFactor, KrylovOrthogonalizeMatrix[k]);
	  }
	NormalizationFactors[i] = KrylovOrthogonalizeMatrix[i] * KrylovOrthogonalizeMatrix[i];
   }
  LongRationalMatrix KrylovInvertOrthogonalizeMatrix = KrylovOrthogonalizeMatrix.DuplicateAndTranspose();
  for (int i = 0; i < KrylovDimension; ++i)
    {
      KrylovOrthogonalizeMatrix[i] /= NormalizationFactors[i];
    }
  // LongRationalMatrix TmpMatrix3 = KrylovOrthogonalizeMatrix * KrylovInvertOrthogonalizeMatrix;
  // if (TmpMatrix3.IsIntegerMatrix())
  //   {
  //     cout << "TmpMatrix3 is an integer matrix" << endl;
  //   }
  cout << "Projecting Hamiltonian" << endl;
  LongRationalMatrix TmpMatrix1 = RationalHamiltonian * KrylovOrthogonalizeMatrix;
  LongRationalMatrix TmpMatrix2 = KrylovInvertOrthogonalizeMatrix * TmpMatrix1;
  
  cout << "Computing characteristic polynomial" << endl;
  LongRational* CharacteristicPolynomial = TmpMatrix2.CharacteristicPolynomial();
  // for (int i = 0; i <= KrylovDimension; ++i)
  //   {
  //     cout << CharacteristicPolynomial[i] << "x^" << i << " + "  << endl;
  //   }

  char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
  sprintf (PolynomialOutputFileName, "%s.charpol", outputFileName);
  ofstream OutputFile;
  OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
  OutputFile << CharacteristicPolynomial[0];
  for (int i = 1; i <= KrylovDimension; ++i)
    {
      OutputFile << "," << CharacteristicPolynomial[i];
    }
  OutputFile << endl;
  OutputFile.close();

  return KrylovDimension;
}

// convert a string describing a product state in the (u,d,+,-) basis to its fermionic occupation version
//
// productState = product state in the (u,d,+,-) basis
// nbrFermions= reference to the number of fermions in the product state
// nbrFermionsEvenMomentum = reference to the sublattice sector
// maxMomentum = reference to the max momentumw number or number of sites
// yMomentum = reference to the y momentumw
// return value = array containing the fermionic occupation (0 if an error occurred) 

int* FQHPairHoppingKrylovParseProductState(char* productState, int& nbrFermions, int& nbrFermionsEvenMomentum, int& maxMomentum, int& yMomentum)
{
  int TmpStringLength = strlen(productState);
  maxMomentum = 2 * TmpStringLength;
  int* ProductStateConfiguration = new int [maxMomentum];
  
  nbrFermions = 0;
  maxMomentum = 0;
  yMomentum = 0;
  nbrFermionsEvenMomentum = 0;
  for (int i = 0; i < TmpStringLength; ++i)
    {
      char TmpChar = productState[i];
      if (TmpChar == '+')
	{
	  ProductStateConfiguration[nbrFermions] = maxMomentum;
	  ProductStateConfiguration[nbrFermions + 1] = maxMomentum + 1;
	  yMomentum += (2 * maxMomentum + 1);
	  nbrFermionsEvenMomentum++;
	  maxMomentum += 2;
	  nbrFermions += 2;
	}
      else
	{
	  if (TmpChar == '-')
	    {
	      maxMomentum += 2;
	    }
	  else
	    {
	      if (TmpChar == 'u')
		{
		  ProductStateConfiguration[nbrFermions] = maxMomentum + 1;
		  yMomentum += (maxMomentum + 1);
		  maxMomentum += 2;
		  nbrFermions++;
		}
	      else
		{
		  if (TmpChar == 'd')
		    {
		      ProductStateConfiguration[nbrFermions] = maxMomentum;
		      nbrFermionsEvenMomentum++;
		      yMomentum += maxMomentum;
		      maxMomentum += 2;
		      nbrFermions++;
		    }
		  else
		    {
		      cout << "illegal character \"" << TmpChar << "\" in " << productState << endl;
		      return 0;
		    }
		}
	    }
	}
    }
  yMomentum %= maxMomentum;
  
  cout << "product state \"" << productState << "\" converted into ";
  for (int i = 0; i < nbrFermions; ++i)
    {
      cout << "c+_{" << ProductStateConfiguration[i] <<"}";
    }
  cout << "|0>" << endl;
  cout << "nbr sites=" << maxMomentum << ", nbr fermions=" << nbrFermions << ", ky=" << yMomentum << ", even sublattice=" << nbrFermionsEvenMomentum << endl;
  return ProductStateConfiguration;
}
