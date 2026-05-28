////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                       DiagHam (FQHE) — N-flavor sphere ED                   //
//                                                                            //
//   New program: FQHESphereFermionsNFlavor                                   //
//   - Supports layer/flavor imbalance via --nbr-particles-per-flavor          //
//   - Reads interaction file keys:                                            //
//       Pseudopotentials_a_b = v0 v1 ... v_{2S}                               //
//       OneBodyPotentials_a_b = t                                       //
//     where a,b run from 1..N (your convention)                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/BooleanOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/StringTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "HilbertSpace/FermionOnSphereWithNFlavor.h"
#include "HilbertSpace/FermionOnSphereWithNFlavorLzSymmetry.h"

#include "Hamiltonian/ParticleOnSphereWithNFlavorGenericHamiltonian.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <stdio.h>

#include <sstream>
#include <vector>
#include <algorithm>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

static void SplitCSVIntegers(const char* s, std::vector<int>& out)
{
  out.clear();
  if (s == 0) return;

  std::string str(s);
  // allow quotes in shell: "2,2,2"
  if (!str.empty() && str.front() == '"') str.erase(str.begin());
  if (!str.empty() && str.back()  == '"') str.pop_back();

  std::stringstream ss(str);
  std::string item;
  while (std::getline(ss, item, ','))
    {
      // trim spaces
      while (!item.empty() && (item.front()==' ' || item.front()=='\t')) item.erase(item.begin());
      while (!item.empty() && (item.back() ==' ' || item.back() =='\t')) item.pop_back();
      if (item.empty()) continue;
      out.push_back(std::atoi(item.c_str()));
    }
}


static inline int ChannelIndex(int s1, int s2, int NbrFlavors)
{
  /*if (NbrFlavors == 2)
  {
    if (s1 == 0 && s2 == 0) return 0;   // up-up
    if (s1 == 1 && s2 == 1) return 1;   // down-down
    return 2;                           // up-down
  }*/
  if (s1 > s2)
    std::swap(s1, s2);

  return s1 * NbrFlavors - (s1*(s1-1))/2 + (s2 - s1);
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager("FQHESphereFermionsNFlavor", "0.01");
  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  

  // -------------------------
  // Option groups
  // -------------------------
  OptionGroup* SystemGroup   = new OptionGroup("system options");
  OptionGroup* ToolsGroup    = new OptionGroup("tools options");
  OptionGroup* OutputGroup   = new OptionGroup("output options");
  OptionGroup* PrecalcGroup  = new OptionGroup("precalculation options");
  OptionGroup* MiscGroup     = new OptionGroup("misc options");

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += PrecalcGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  // -------------------------
  // System options (new ones)
  // -------------------------
  (*SystemGroup) += new SingleIntegerOption('f', "nbr-flavors", "number of flavors (layers)", 2);
  (*SystemGroup) += new SingleStringOption('\n', "nbr-particles-per-flavor",
                                           "comma-separated list n0,n1,...,n(NbrFlavors-1) (must have NbrFlavors entries)",
                                           0);

  // Keep common/legacy options (so you can reuse your existing workflow)
  (*SystemGroup) += new SingleIntegerOption('p', "nbr-particles", "total number of particles (optional sanity check)", 0);
  (*SystemGroup) += new SingleIntegerOption('l', "lzmax", "twice the maximum single particle angular momentum (2S)", 8);
  (*SystemGroup) += new SingleIntegerOption('\n', "total-lz", "twice total Lz (sector)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);

  // Lz <-> -Lz symmetry sector machinery (only valid if total-lz = 0)
  (*SystemGroup) += new BooleanOption('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized basis (only if total-lz=0)", false);
  (*SystemGroup) += new BooleanOption('\n', "minus-lzparity", "select negative parity sector for Lz symmetry (only if lzsymmetrized-basis)", false);

  // Interaction file options
  (*SystemGroup) += new SingleStringOption('\n', "interaction-file", "interaction (pseudopotential) file name", 0);
  (*SystemGroup) += new SingleStringOption('\n', "interaction-name", "interaction name label used in output files", (char*)"nflavor");
  (*SystemGroup) += new SingleStringOption('\n', "use-hilbert", "external Hilbert space file", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "fast-search", "memory for Hilbert space lookup (MB)", 100);

  (*OutputGroup) += new BooleanOption('\n', "eigenstate", "compute and store eigenvectors", false);
  (*OutputGroup) += new BooleanOption('\n', "use-entanglement", "enable entanglement-related outputs (if available)", false);
  (*OutputGroup) += new SingleDoubleOption('\n', "energy-shift", "energy shift applied to Hamiltonian (for Lanczos stability)", 0.0);

  (*PrecalcGroup) += new BooleanOption('\n', "disk-cache", "use on-disk cache for precalculations", false);
  (*PrecalcGroup) += new SingleIntegerOption('m', "memory", "precalculation memory in MB", 0);
  (*PrecalcGroup) += new SingleStringOption('\n', "save-precalculation", "save precalculation to this file", 0);
  (*PrecalcGroup) += new SingleStringOption('\n', "load-precalculation", "load precalculation from this file", 0);
  
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption('\n', "use-lapack", "use LAPACK instead of DiagHam");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption('h', "help", "display this help", false);

  // -------------------------
  // Parse options
  // -------------------------
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp(cout);
      return 0;
    }

  const int NbrFlavors= Manager.GetInteger("nbr-flavors");
  const int LzMax = Manager.GetInteger("lzmax");
  const int TotalLz = Manager.GetInteger("total-lz");
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();

  bool FirstRun = true;

  const char* NPerFlavorStr = Manager.GetString("nbr-particles-per-flavor");
  if (NPerFlavorStr == 0)
    {
      cerr << "Error: missing --nbr-particles-per-flavor \"n0,n1,...\" " << endl;
      return -1;
    }

  std::vector<int> NPerFlavor;
  SplitCSVIntegers(NPerFlavorStr, NPerFlavor);

  if ((int)NPerFlavor.size() != NbrFlavors )
    {
      cerr << "Error: --nbr-particles-per-flavor must contain exactly NbrFlavors =" << NbrFlavors 
           << " integers, got " << (int)NPerFlavor.size() << endl;
      return -1;
    }

  int NbrFermions = 0;
  for (int i=0; i<NbrFlavors; ++i)
    {
      if (NPerFlavor[i] < 0)
        {
          cerr << "Error: negative particle number in --nbr-particles-per-flavor." << endl;
          return -1;
        }
      NbrFermions += NPerFlavor[i];
    }

  // optional sanity check
  const int NbrFermionsOpt = Manager.GetInteger("nbr-particles");
  if (NbrFermionsOpt > 0 && NbrFermionsOpt != NbrFermions)
    {
      cerr << "Error: --nbr-particles (" << NbrFermionsOpt << ") != sum(--nbr-particles-per-flavor) ("
           << NbrFermions << ")" << endl;
      return -1;
    }

  cout << "Populations: ";
  for (int i = 0; i < NbrFlavors ; ++i)
  cout << "N" << i << "=" << NPerFlavor[i] << " ";
  cout << endl;
  
  const bool UseLzSym = Manager.GetBoolean("lzsymmetrized-basis");
  const bool MinusParity = Manager.GetBoolean("minus-lzparity");
  if (UseLzSym && TotalLz != 0)
    {
      cerr << "Error: --lzsymmetrized-basis is only valid for --total-lz 0" << endl;
      return -1;
    }

  const char* InteractionFileName = Manager.GetString("interaction-file");
  if (InteractionFileName == 0)
    {
      cerr << "Error: missing --interaction-file" << endl;
      return -1;
    }

  const char* InteractionName = Manager.GetString("interaction-name");
  if (InteractionName == 0) InteractionName = "nflavor";

  // -------------------------
  // Read interaction file: Pseudopotentials_a_b and OneBodyPotentials_a_b
  // -------------------------
  ConfigurationParser InteractionFile;
  if (InteractionFile.Parse(InteractionFileName) == false)
    {
      cerr << "Error: cannot read interaction file '" << InteractionFileName << "'" << endl;
      return -1;
    }

  // Two-body: store fullNbrFlavorsxNbrFlavors(we’ll just fill what’s provided; missing -> 0)
  int NChannels = NbrFlavors *(NbrFlavors+1)/2;
  double** TwoBody = new double*[NChannels];
  for (int ab=0; ab<NChannels; ++ab)
    {
      TwoBody[ab] = new double[LzMax+1];
      for (int m=0;m<=LzMax;++m) TwoBody[ab][m] = 0.0;
    }

  // One-body tunneling: full NbrFlavors x NbrFlavors single number (missing -> 0)
  double*** OneBody = new double**[NbrFlavors];
  for (int a=0; a<NbrFlavors; ++a)
  {
      OneBody[a] = new double*[NbrFlavors];
      for (int b=0; b<NbrFlavors; ++b)
      {
          OneBody[a][b] = new double[LzMax+1];
          for (int m=0; m<=LzMax; ++m)
              OneBody[a][b][m] = 0.0;
      }
  }

  // Parse all Pseudopotentials_a_b
  for (int a=1; a<=NbrFlavors; ++a)
    for (int b=1; b<=NbrFlavors; ++b)
      {
        char key[256];
        int channel = ChannelIndex(a-1,b-1, NbrFlavors);  // 0-based

        // Two-body
        std::sprintf(key, "Pseudopotentials_%d_%d", a, b);
        double* vals = 0;
        int nvals = 0;
        if (InteractionFile.GetAsDoubleArray(key, ' ', vals, nvals) == true)
          {
            int copyMax = (nvals <= (LzMax+1)) ? nvals : (LzMax+1);
            for (int m = 0; m < copyMax; ++m)
              TwoBody[channel][m] = vals[m];
            delete[] vals;
          }

        // One-body
        std::sprintf(key, "OneBodyPotentials_%d_%d", a, b);
        vals = 0;
        nvals = 0;
        if (InteractionFile.GetAsDoubleArray(key, ' ', vals, nvals) == true)
          {
            int copyMax = (nvals <= (LzMax+1)) ? nvals : (LzMax+1);
            for (int m=0; m<copyMax; ++m)
              OneBody[a-1][b-1][m] = vals[m];
            delete[] vals;
          }
      }


  //char OutputNameLz[512];
  char* OutputNameLz = new char [512 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf(OutputNameLz, "fermions_sphere_suN%d_%s_n_%d_2s_%d_lz.dat", NbrFlavors, Manager.GetString("interaction-name"), NbrFermions, LzMax);

  // -------------------------
  // Build Hilbert space (with imbalance) and Hamiltonian
  // -------------------------
  
  unsigned long Memory = ((unsigned long)Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long)Manager.GetInteger("fast-search")) << 20;
  
  const bool onDiskCacheFlag = Manager.GetBoolean("disk-cache");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");

  
  int Max = 0;
  for (int f = 0; f < NbrFlavors; ++f)
    Max += ((LzMax - NPerFlavor[f] + 1) * NPerFlavor[f]);
  cout << "maximum Lz value = " << Max << endl;

  int L = 0;

  if ((abs(Max) & 1) != 0)
    L = 1;

  if (InitialLz >= 0)
  {
  L = InitialLz;
  if ((abs(Max) & 1) != 0)
    L |= 1;
  else
    L &= ~0x1;
  }

  if (NbrLz > 0)
  {
  if (L + (2 * (NbrLz - 1)) < Max)
    Max = L + (2 * (NbrLz - 1));
  }
  
  cout << "LzMax = " << LzMax << endl;
  cout << "Total flavors = " << NbrFlavors << endl;
  cout << "Total modes = " << (LzMax+1)*NbrFlavors << endl;
  cout << "Lz sectors from " << L << " to " << Max << endl;

  //ParticleOnSphere* Space = 0;
  //if (UseLzSym)
  //  { Space = new FermionOnSphereWithNFlavorLzSymmetry( NbrFermions, TotalLz, LzMax, N, &(NPerFlavor), Memory
                //N, &(NPerFlavor[0]), LzMax,
                //MinusParity,
                //Memory
    //          );
    //}
  //else
  //  { Space = new FermionOnSphereWithNFlavor( LzMax, TotalLz, N, NPerFlavor
                //N, &(NPerFlavor[0]), LzMax,
                //TotalLz,
                //Memory );
  //  }

 for (; L <= Max; L += 2)
 {
  //cout << "Running Lz sector = " << L << endl;

  FermionOnSphereWithNFlavor* Space = 0;
  if (UseLzSym)
    Space = new FermionOnSphereWithNFlavorLzSymmetry(NbrFermions, L, LzMax, NbrFlavors, &(NPerFlavor[0]), MemorySpace);
  else
    Space = new FermionOnSphereWithNFlavor(NbrFermions, L, LzMax, NbrFlavors, &(NPerFlavor[0]), MemorySpace);

  //cout << "Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

  //Architecture.SetDimension(Space->GetHilbertSpaceDimension());
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();

    //cout << "Program reached Hilbert Space safely." << endl;

  AbstractQHEHamiltonian* Hamiltonian = 0;
  Hamiltonian = new ParticleOnSphereWithNFlavorGenericHamiltonian(Space, NbrFermions, LzMax, TwoBody, OneBody, Architecture.GetArchitecture(), Memory, onDiskCacheFlag, SavePrecalculationFileName);

  //Hamiltonian = new ParticleOnSphereWithNFlavorGenericHamiltonian( (FermionOnSphereWithNFlavor*) Space, LzMax, N, InteractionFileName );
  //Hamiltonian = new ParticleOnSphereWithNFlavorGenericHamiltonian(Space, LzMax, N, InteractionFileName);

  //const double Shift = Manager.GetDouble("energy-shift");
  double Shift = -10.0;
  if (Shift != 0.0) Hamiltonian->ShiftHamiltonian(Shift);

  if (SavePrecalculationFileName != 0)
      Hamiltonian->SavePrecalculation(SavePrecalculationFileName);

    //cout << "Program reached Hamiltonian safely." << endl;
  // -------------------------
  // Run the standard sphere main task (Lanczos, outputs, etc.)
  // -------------------------
  
    char* EigenvectorName = 0;
    if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
    {
      EigenvectorName = new char [120];
      sprintf (EigenvectorName, "fermions_sphere_suN%d_%s_n_%d_2s_%d_lz_%d", NbrFlavors, ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, LzMax, L);
    }
    QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
    MainTaskOperation TaskOperation (&Task);
    TaskOperation.ApplyOperation(Architecture.GetArchitecture());
    //Task.ExecuteMainTask();

    delete Hamiltonian;
    delete Space;
    if (EigenvectorName != 0)
	  {
	    delete[] EigenvectorName;
	    EigenvectorName = 0;
	  }

    if (FirstRun == true) FirstRun = false;
  }
  // -------------------------
  // Cleanup
  // -------------------------


  for (int ab=0; ab<NChannels; ++ab) delete[] TwoBody[ab];
  delete[] TwoBody;
  for (int a=0; a<NbrFlavors; ++a)
  {
    for (int b=0; b<NbrFlavors; ++b)
      delete[] OneBody[a][b];
    delete[] OneBody[a];
  }
  delete[] OneBody;

  delete[] OutputNameLz;
  return 0;
}
