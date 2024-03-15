#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#ifdef __GSL__
#include <gsl/gsl_sf_gamma.h>
#endif

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital
// nbrFluxQuanta  = number of flux quanta
// theta = polar angle that defines the cut (in pi units)
// return value = square of the orbital weight
double FQHESphereComputeSharpRealSpaceCutCoefficient (int orbitalIndex, int nbrFluxQuanta, double theta);

// compute the weight of a given orbital for a sharp real space cut for a finite patch both in theta and phi
//
// orbitalIndex1 = index of the first orbital
// orbitalIndex2 = index of the second orbital
// nbrFluxQuanta  = number of flux quanta
// theta = polar angle that defines the cut (in pi units)
// phi = azimuthal angle that defines the cut (in pi units)
// return value = square of the orbital weight
double FQHESphereComputeSharpRealSpaceCutCoefficient (int orbitalIndex1, int orbitalIndex2, int nbrFluxQuanta, double theta, double phi);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereRealSpacePartitionCoefficients" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta", 10);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "theta", "polar angle that defines the cut (in pi units)", 0.5);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "phi", "azimuthal angle that defines the cut (in pi units). The region is defined between -phi/2 and phi/2. 0 preserves the rotation symmetry along z (i.e. equivalenet to 2pi)", 0.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "optional output file name (default is realspace_disk_theta_*_phi_*_2s_*.dat)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereRealSpacePartitionCoefficients -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");

  double Theta = Manager.GetDouble("theta");
  double Phi = Manager.GetDouble("phi");
  double CutPosition = 0.0;
  char* OutputFile = 0;
  if (Manager.GetString("output-file") == 0)
    {
      OutputFile = new char[512];
      if (Phi == 0.0)
	{
	  sprintf (OutputFile, "realspace_sphere_theta_%.6f_2s_%d.dat", Theta, NbrFluxQuanta);
	}
      else
	{
	  sprintf (OutputFile, "realspace_sphere_theta_%.6f_phi_%.6f_2s_%d.dat", Theta, Phi, NbrFluxQuanta);
	}
    }
  else
    {
      OutputFile = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFile, Manager.GetString("output-file"));
    }
  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  if (Phi == 0.0)
    {
      File << "# real space coefficients for a cut at theta=" << Theta << " on a sphere with N_phi=" << NbrFluxQuanta << endl
	   << "OrbitalSquareWeights =";
      int NbrCoefficients = 0;
      double* Coefficients = 0;
      Coefficients = new double [NbrFluxQuanta + 1];
      for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
	{
	  Coefficients[NbrCoefficients] = FQHESphereComputeSharpRealSpaceCutCoefficient(NbrCoefficients, NbrFluxQuanta, Theta);
	}
      for (int i = 0; i < NbrCoefficients; ++i)
	File << " " << Coefficients[i];
      File << endl;
      delete[] Coefficients;
    }
  else
    {
      File << "# real space coefficients for a cut at theta=" << Theta << " and phi=" << Phi << " on a sphere with N_phi=" << NbrFluxQuanta << endl;
      RealSymmetricMatrix TmpOverlapMatrix (NbrFluxQuanta + 1, true);
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  TmpOverlapMatrix.SetMatrixElement(i, i, FQHESphereComputeSharpRealSpaceCutCoefficient(i, i, NbrFluxQuanta, Theta, Phi));
	  for (int j = i + 1; j <= NbrFluxQuanta; ++j)
	    {
	      TmpOverlapMatrix.SetMatrixElement(i, j, FQHESphereComputeSharpRealSpaceCutCoefficient(i, j, NbrFluxQuanta, Theta, Phi));
	    }
	}
      
      RealMatrix TmpTransformationMatrix1 (NbrFluxQuanta + 1, NbrFluxQuanta + 1);
      TmpTransformationMatrix1.SetToIdentity();
      RealDiagonalMatrix TmpDiag1 = RealDiagonalMatrix(NbrFluxQuanta + 1, true);  
      RealDiagonalMatrix TmpDiag2 = RealDiagonalMatrix(NbrFluxQuanta + 1, true);  
#ifdef __LAPACK__
      TmpOverlapMatrix.LapackDiagonalize(TmpDiag1, TmpTransformationMatrix1);
#else
      TmpOverlapMatrix.Diagonalize(TmpDiag1, TmpTransformationMatrix1);
#endif 	  
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  //	  cout << TmpDiag1[i] << endl;
	  TmpDiag2[i] = sqrt(1.0 - TmpDiag1[i]);	      
	  TmpDiag1[i] = sqrt(TmpDiag1[i]);
	}
      RealMatrix TmpTransformationMatrix2 = TmpTransformationMatrix1.DuplicateAndTranspose();
      RealMatrix TmpTransformationMatrix3 = TmpTransformationMatrix1 * TmpDiag2;
      for (int i = 0; i < TmpTransformationMatrix1.GetNbrColumn(); ++i)
	{
	  TmpTransformationMatrix1[i] *= TmpDiag1[i];
	}
      RealMatrix TmpTransformationMatrix4 = TmpTransformationMatrix1 * TmpTransformationMatrix2;
      RealMatrix TmpTransformationMatrix5 = TmpTransformationMatrix3 * TmpTransformationMatrix2;
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  for (int j = 0; j <= NbrFluxQuanta; ++j)
	    {
	      double Tmp1;
	      double Tmp2;
	      TmpTransformationMatrix4.GetMatrixElement(i, j, Tmp1);
	      TmpTransformationMatrix5.GetMatrixElement(i, j, Tmp2);
	      File << i << " " << j << " " << Tmp1 << " " << Tmp2 << endl;
	    }
	}
    }
  File.close();
  return 0;
}


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital
// nbrFluxQuanta  = number of flux quanta
// theta = polar angle that defines the cut (int pi units)
// return value = square of the orbital weight

double FQHESphereComputeSharpRealSpaceCutCoefficient (int orbitalIndex, int nbrFluxQuanta, double theta)
{  
#ifdef __GSL__
  return gsl_sf_beta_inc((double) (orbitalIndex + 1), (double) (nbrFluxQuanta - orbitalIndex + 1), sin(theta * M_PI * 0.5) * sin(theta * M_PI * 0.5));
#else
  cout << "gsl is required" << endl;
  return 0.0;
#endif
}

// compute the weight of a given orbital for a sharp real space cut for a finite patch both in theta and phi
//
// orbitalIndex1 = index of the first orbital
// orbitalIndex2 = index of the second orbital
// nbrFluxQuanta  = number of flux quanta
// theta = polar angle that defines the cut (in pi units)
// phi = azimuthal angle that defines the cut (in pi units)
// return value = square of the orbital weight

double FQHESphereComputeSharpRealSpaceCutCoefficient (int orbitalIndex1, int orbitalIndex2, int nbrFluxQuanta, double theta, double phi)
{
#ifdef __GSL__
  if (orbitalIndex1 == orbitalIndex2)
    {
      return (gsl_sf_beta_inc((double) (orbitalIndex1 + 1), (double) (nbrFluxQuanta - orbitalIndex1 + 1), sin(theta * M_PI * 0.5) * sin(theta * M_PI * 0.5)) * phi);
    }
  else
    {
      return ((gsl_sf_beta_inc((0.5 * ((double) (orbitalIndex1 + orbitalIndex2))) + 1.0,
			      (0.5 * ((double) ((2 * nbrFluxQuanta) - (orbitalIndex1 + orbitalIndex2)))) + 1.0,
			      sin(theta * M_PI * 0.5) * sin(theta * M_PI * 0.5)) * (sin(((double) (orbitalIndex1 - orbitalIndex2)) * M_PI * phi)) / (M_PI * ((double) (orbitalIndex1 - orbitalIndex2))))
	      * (gsl_sf_beta((0.5 * ((double) (orbitalIndex1 + orbitalIndex2))) + 1.0,
			     (0.5 * ((double) ((2 * nbrFluxQuanta) - (orbitalIndex1 + orbitalIndex2)))) + 1.0) / sqrt(gsl_sf_beta((double) (orbitalIndex1 + 1), (double) (nbrFluxQuanta - orbitalIndex1 + 1)) * gsl_sf_beta((double) (orbitalIndex2 + 1), (double) (nbrFluxQuanta - orbitalIndex2 + 1)))));
    }
#else
  cout << "gsl is required" << endl;
  return 0.0;
#endif
}

