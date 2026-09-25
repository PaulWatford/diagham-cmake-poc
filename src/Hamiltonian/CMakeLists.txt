# Auto-generated CMakeLists.txt for src/Hamiltonian
# Regenerate with: python3 scripts_cmake/extract_autotools.py

add_subdirectory(DMRGHamiltonian)

diagham_add_library(Hamiltonian
    SOURCES
        AbstractHamiltonian.cc
        ExplicitHamiltonian.cc
        FileBasedHamiltonian.cc
        FileBasedHermitianHamiltonian.cc
        TensorProductSparseMatrixHamiltonian.cc
        TensorProductSparseMatrixSelectedBlockHamiltonian.cc
        MultipleTensorProductSparseMatrixHamiltonian.cc
        TripleTensorProductSparseMatrixHamiltonian.cc
)
