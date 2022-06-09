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


// compute the characteristic polynomial for hamiltonian projected onto the Krylov subspace
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// productStateIndex = index of the product state that defines the Krylov subpspace
// return value = Krylov subspace dimension
int FQHPairHoppingKrylovSubspaceComputeCharacteristic(AbstractQHEHamiltonian* hamiltonian, ParticleOnTorusWithMagneticTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, int productStateIndex);


int main(int argc, char** argv)
{
  cout.precision(14);

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
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
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

  if (Manager.GetString("product-state") == 0)
    {
      cout << "error, a product state should be provided" << endl;
      return 0;
    }

  int TmpStringLength = strlen(Manager.GetString("product-state"));
  MaxMomentum = 2 * TmpStringLength;
  int* ProductStateConfiguration = new int [MaxMomentum];
  
  NbrFermions = 0;
  MaxMomentum = 0;
  YMomentum = 0;
  NbrFermionsEvenMomentum = 0;
  for (int i = 0; i < TmpStringLength; ++i)
    {
      char TmpChar = Manager.GetString("product-state")[i];
      if (TmpChar == '+')
	{
	  ProductStateConfiguration[NbrFermions] = MaxMomentum;
	  ProductStateConfiguration[NbrFermions + 1] = MaxMomentum + 1;
	  YMomentum += (2 * MaxMomentum + 1);
	  NbrFermionsEvenMomentum++;
	  MaxMomentum += 2;
	  NbrFermions += 2;
	}
      else
	{
	  if (TmpChar == '-')
	    {
	      MaxMomentum += 2;
	    }
	  else
	    {
	      if (TmpChar == 'u')
		{
		  ProductStateConfiguration[NbrFermions] = MaxMomentum + 1;
		  YMomentum += (MaxMomentum + 1);
		  MaxMomentum += 2;
		  NbrFermions++;
		}
	      else
		{
		  if (TmpChar == 'd')
		    {
		      ProductStateConfiguration[NbrFermions] = MaxMomentum;
		      NbrFermionsEvenMomentum++;
		      YMomentum += MaxMomentum;
		      MaxMomentum += 2;
		      NbrFermions++;
		    }
		  else
		    {
		      cout << "illegal character \"" << TmpChar << "\" in " << Manager.GetString("product-state") << endl;
		      return 0;
		    }
		}
	    }
	}
    }
  YMomentum %= MaxMomentum;
  
  cout << "product state \"" << Manager.GetString("product-state") << "\" converted into ";
  for (int i = 0; i < NbrFermions; ++i)
    {
      cout << "c+_{" << ProductStateConfiguration[i] <<"}";
    }
  cout << "|0>" << endl;
  cout << "nbr sites=" << MaxMomentum << ", nbr fermions=" << NbrFermions << ", ky=" << YMomentum << ", even sublattice=" << NbrFermionsEvenMomentum << endl;
  
  
  char* SuffixOutputName = new char [256];
  sprintf (SuffixOutputName, "n_%d_2s_%d.dat", NbrFermions, MaxMomentum);

  char* OutputName = new char [512 + strlen(SuffixOutputName)];
  sprintf (OutputName, "fermions_pairhopping_krylov_%s_%s", Manager.GetString("product-state"), SuffixOutputName);

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
	  int ProductStateIndex = TotalSpace->FindStateIndex(ProductStateConfiguration);
	  if (ProductStateIndex == -1)
	    {
	      cout << "error, the product state configuration \"" << Manager.GetString("product-state") << "\" is not compatible with the Hilbert space" << endl;
	      return 0;
	    }
	  cout << "Product state configuration converted as ";
	  TotalSpace->PrintState(cout, ProductStateIndex) << " (index=" << ProductStateIndex << ")" << endl;
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
	  
	  
	  int KrylovSpaceDimension = FQHPairHoppingKrylovSubspaceComputeCharacteristic(Hamiltonian, TotalSpace, EigenvectorName, Architecture.GetArchitecture(), ProductStateIndex);

	  if (Manager.GetBoolean("use-double"))
	    {
	      RealVector* TmpVectors = new RealVector[TotalSpace->GetHilbertSpaceDimension()];
	      TmpVectors[0] = RealVector(TotalSpace->GetHilbertSpaceDimension(), true);
	      TmpVectors[0][ProductStateIndex] = 1.0;
	      FullReorthogonalizedLanczosAlgorithm TmpLanczos (Architecture.GetArchitecture(), KrylovSpaceDimension, KrylovSpaceDimension + 1);
	      TmpLanczos.SetHamiltonian(Hamiltonian);
	      TmpLanczos.InitializeLanczosAlgorithm(TmpVectors[0]);
	      TmpLanczos.RunLanczosAlgorithm(KrylovSpaceDimension);
	      RealVector* KrylovLanczosVectors = (RealVector*) TmpLanczos.GetKrylovSubspace();
	      RealTriDiagonalSymmetricMatrix TmpTridiagHamiltonian = TmpLanczos.GetTridiagonalMatrix();
	      RealMatrix TridiagHamiltonian (TmpTridiagHamiltonian);
	      double* TridiagCharacteristicPolynomial = TridiagHamiltonian.CharacteristicPolynomial();
	      for (int i = 0; i <= KrylovSpaceDimension; ++i)
		{
		  cout << TridiagCharacteristicPolynomial[i] << "x^" << i << " + "  << endl;
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
// productStateIndex = index of the product state that defines the Krylov subpspace
// return value = Krylov subspace dimension

int FQHPairHoppingKrylovSubspaceComputeCharacteristic(AbstractQHEHamiltonian* hamiltonian, ParticleOnTorusWithMagneticTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, int productStateIndex)
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
#ifdef __GMP__
  mpz_set_ui(TmpVectors[0][productStateIndex], 1ul);
#else
  TmpVectors[0][productStateIndex] = 1l;  
#endif
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
