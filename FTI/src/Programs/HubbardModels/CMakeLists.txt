# Auto-generated CMakeLists.txt for FTI/src/Programs/HubbardModels
# Regenerate with: python3 scripts_cmake/extract_autotools.py

# Programs listed in the upstream Makefile.am bin_PROGRAMS.
diagham_add_programs(PREFIX HubbardModels
    PROGRAMS
        HubbardShowBasis HubbardShowBasis.cc
        HubbardGetDimension HubbardGetDimension.cc
        HubbardEntanglementEntropyParticlePartition HubbardEntanglementEntropyParticlePartition.cc
        HubbardEntanglementEntropyRealSpacePartition HubbardEntanglementEntropyRealSpacePartition.cc
        HubbardKitaevHeisenbergModel HubbardKitaevHeisenbergModel.cc
        HubbardSuperconductorOrderParameter HubbardSuperconductorOrderParameter.cc
        HubbardASCIIToState HubbardASCIIToState.cc
        HubbardConvertFromTranslationInvariantBasis HubbardConvertFromTranslationInvariantBasis.cc
        HubbardDensity HubbardDensity.cc
        HubbardKitaevHeisenbergHoneycombModel HubbardKitaevHeisenbergHoneycombModel.cc
        HubbardFloquetHeisenbergModel HubbardFloquetHeisenbergModel.cc
        HubbardSquareLatticeModel HubbardSquareLatticeModel.cc
        HubbardSquareLatticeEtaPairingStateTimeEvolution HubbardSquareLatticeEtaPairingStateTimeEvolution.cc
        HubbardKagomeLatticeModel HubbardKagomeLatticeModel.cc
        HubbardSquareLatticeEtaPairingEntanglementEntropyRealSpacePartition HubbardSquareLatticeEtaPairingEntanglementEntropyRealSpacePartition.cc
        HubbardSquareLatticeModelJ2S2 HubbardSquareLatticeModelJ2S2.cc
        HubbardApplyGutzwillerProjection HubbardApplyGutzwillerProjection.cc
        KitaevChain KitaevChain.cc
        MultipleKitaevChainsToSpin1 MultipleKitaevChainsToSpin1.cc
        HubbardHaldaneLatticeModel HubbardHaldaneLatticeModel.cc
        HubbardSSHModel HubbardSSHModel.cc
        HubbardSSHModelMPS HubbardSSHModelMPS.cc
        HubbardExtendedKitaevHeisenbergHoneycombModel HubbardExtendedKitaevHeisenbergHoneycombModel.cc
        HubbardHoneycombClusterChargingLatticeModel HubbardHoneycombClusterChargingLatticeModel.cc
        KondoModelFlatBands KondoModelFlatBands.cc
)