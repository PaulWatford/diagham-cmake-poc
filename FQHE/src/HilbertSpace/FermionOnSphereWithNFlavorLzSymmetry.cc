#include "config.h"
#include "HilbertSpace/FermionOnSphereWithNFlavorLzSymmetry.h"

FermionOnSphereWithNFlavorLzSymmetry::FermionOnSphereWithNFlavorLzSymmetry(
    int nbrFermions,
    int totalLz,
    int lzMax,
    int nbrFlavors,
    int* nbrParticlesPerFlavor,
    unsigned long memory)
  : FermionOnSphereWithNFlavor(nbrFermions, totalLz, lzMax, nbrFlavors, nbrParticlesPerFlavor, memory)
{
}

FermionOnSphereWithNFlavorLzSymmetry::~FermionOnSphereWithNFlavorLzSymmetry()
{
}

int FermionOnSphereWithNFlavorLzSymmetry::GetHilbertSpaceAdditionalSymmetry()
{
  // 0 = "no extra symmetry sector selected" (placeholder)
  // You can later encode +/- Lz parity here like SpinLzSymmetry does.
  return 0;
}


AbstractHilbertSpace*
FermionOnSphereWithNFlavorLzSymmetry::Clone()
{
    return new FermionOnSphereWithNFlavorLzSymmetry(*this);
}

std::ostream&
FermionOnSphereWithNFlavorLzSymmetry::PrintState(std::ostream& os, int state)
{
    return FermionOnSphereWithNFlavor::PrintState(os, state);
}

int
FermionOnSphereWithNFlavorLzSymmetry::GetParticleStatistic()
{
    return FermionicStatistic;
}