#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong.h"

#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/FilenameTools.h"

#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FTIGenerateMagnon" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  ArchitectureManager Architecture;
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "name of the vector file to which the magnon operator should be applied");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-groundstate", "name of the file that gives the vector files to which the  magnon operator should be applied (in all-bilinear or sma mode)");	
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  //  (*SystemGroup) += new BooleanOption  ('\n', "all-bilinear", "apply all bilinear operators to the ground state, without summing on the sector");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIGenerateMagnon -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if ((Manager.GetString("eigenstate-file") == 0) && (Manager.GetString("degenerate-groundstate") == 0))
    {
      cout << "error, an eigenstate state file should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
  
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int* SzValues = 0;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  bool Statistics = true;
  ComplexVector* GroundStates;
  ParticleOnSphereWithSpin** SpaceSource = 0;

  if (Manager.GetString("degenerate-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      SzValues = new int[1];
      GroundStates = new ComplexVector[1];
      SpaceSource = new ParticleOnSphereWithSpin* [1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("eigenstate-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("eigenstate-file"));
      if (GroundStates[0].ReadVector (GroundStateFiles[0]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[0] << endl;
	  return -1;      
	}	
      if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSitesX, NbrSitesY, TotalKx[0], TotalKy[0], SzValues[0], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
	  return -1;
	}
      cout << GroundStateFiles[0] << " has quantum numbers N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " Kx=" << TotalKx[0] << " Ky=" << TotalKy[0] << " 2Sz=" << SzValues[0] << endl;
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerate-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
      NbrSpaces = DegeneratedFile.GetNbrLines();
      GroundStateFiles = new char* [NbrSpaces];
      TotalKx = new int[NbrSpaces];
      TotalKy = new int[NbrSpaces];
      SzValues = new int[NbrSpaces];
      GroundStates = new ComplexVector[NbrSpaces];
      SpaceSource = new ParticleOnSphereWithSpin* [NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));	
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }	
	  if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i], SzValues[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	  cout << GroundStateFiles[i] << " has quantum numbers N=" << NbrParticles << " Nx=" << NbrSitesX << " Ny=" << NbrSitesY << " Kx=" << TotalKx[i] << " Ky=" << TotalKy[i] << " 2Sz=" << SzValues[i] << endl;
	}
    }
  
  if (Statistics == true)
    {
      if ((NbrSitesX * NbrSitesY) <= 32)
	{
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      SpaceSource[i] = new FermionOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (NbrParticles + SzValues[i]) >> 1, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
	    }
	}
      else
	{
    	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      SpaceSource[i] = new FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong (NbrParticles, (NbrParticles + SzValues[i]) >> 1, NbrSitesX, NbrSitesY, TotalKx[i], TotalKy[i]);
	    }
	}
    }
  else
    {
      cout << "bosons are not implemented yet" << endl;
      return 0;
    }

  
  double* DummyChemicalPotentials = new double[1];
  DummyChemicalPotentials[0] = 0.0;
  Abstract2DTightBindingModel* TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, 1, DummyChemicalPotentials,
											      0.0, 0.0, Architecture.GetArchitecture(), true);
  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
 
  for (; MinKx <= MaxKx; ++MinKx)
    {
      int MinKy = 0;
      int MaxKy = NbrSitesY - 1;
      if (Manager.GetInteger("only-ky") >= 0)
	{						
	  MinKy = Manager.GetInteger("only-ky");
	  MaxKy = MinKy;
	}
      for (; MinKy <= MaxKy; ++MinKy)
	{
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      ParticleOnSphereWithSpin* SpaceDestination = 0;
	      int TargetSz = SzValues[i];
	      int AnnihilationSpinIndex = 0;
	      int CreationSpinIndex = 1;
	      if (TargetSz > 0)
		{
		  TargetSz -= 2;
		  AnnihilationSpinIndex = 1;
		  CreationSpinIndex = 0;
		}
	      else
		{
		  TargetSz += 2;
		}
	      if (Statistics == true)
		{
		  if ((NbrSitesX * NbrSitesY) <= 32)
		    {
		      SpaceDestination = new FermionOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (NbrParticles + TargetSz) >> 1, NbrSitesX, NbrSitesY, MinKx, MinKy);
		    }
		  else
		    {
		      SpaceDestination = new FermionOnSquareLatticeWithSU2SpinMomentumSpaceLong (NbrParticles, (NbrParticles + TargetSz) >> 1, NbrSitesX, NbrSitesY, MinKx, MinKy);
		    }
		}
	      SpaceSource[i]->SetTargetSpace(SpaceDestination);
	      cout << "From Hilbert space dim=" << SpaceSource[i]->GetHilbertSpaceDimension() << " to Hilbert space dim=" << SpaceDestination->GetHilbertSpaceDimension() << endl;
	      for (int Q1x = 0; Q1x < NbrSitesX; ++Q1x)
		{
		  for (int Q1y = 0; Q1y < NbrSitesY; ++Q1y)
		    {
		      int Q2x = TotalKx[i] + Q1x - MinKx;
		      while (Q2x < 0)
			{
			  Q2x += NbrSitesX;
			}
		      Q2x %= NbrSitesX;
		      int Q2y = TotalKy[i] + Q1y - MinKy;
		      while (Q2y < 0)
			{
			  Q2y += NbrSitesY;
			}
		      Q2y %= NbrSitesY;
		      ComplexVector EigenstateOutput(SpaceDestination->GetHilbertSpaceDimension(), true);
		      ParticleOnSphereWithSpinDensityOperator Projector(SpaceSource[i], TightBindingModel->GetLinearizedMomentumIndex(Q1x, Q1y), CreationSpinIndex,
									TightBindingModel->GetLinearizedMomentumIndex(Q2x, Q2y), AnnihilationSpinIndex);
		      VectorOperatorMultiplyOperation Operation(&Projector, &(GroundStates[i]), &EigenstateOutput);
		      Operation.ApplyOperation(Architecture.GetArchitecture());
		      double TmpNorm = EigenstateOutput.Norm();
		      cout << "norm = " << TmpNorm << endl;
		      if (TmpNorm != 0.0)
			{
			  EigenstateOutput /= TmpNorm;
			  char* TmpExtention = new char[128];
			  sprintf (TmpExtention, "_magnon_qx_%d_qy_%d_q1x_%d_q1y_%d_bz_%d.vec", MinKx, MinKy, Q2x, Q2y, TargetSz);
			  char* EigenstateOutputFile = ReplaceExtensionToFileName(GroundStateFiles[i], ".vec", TmpExtention);
			  EigenstateOutput.WriteVector(EigenstateOutputFile);
			  delete[] TmpExtention;
			  delete[] EigenstateOutputFile;
			}
		    }
		}
	    }
	}
    }
  return 0;
}
