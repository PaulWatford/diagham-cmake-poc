# Auto-generated CMakeLists.txt for src/Matrix
# Regenerate with: python3 scripts_cmake/extract_autotools.py

diagham_add_library(Matrix
    SOURCES
        BlockDiagonalMatrix.cc
        BlockMatrix.cc
        ComplexMatrix.cc
        ComplexTriDiagonalHermitianMatrix.cc
        ComplexSkewSymmetricMatrix.cc
        HermitianMatrix.cc
        Matrix.cc
        RealAntisymmetricMatrix.cc
        RealDiagonalMatrix.cc
        RealMatrix.cc
        RealSymmetricMatrix.cc
        RealTriDiagonalSymmetricMatrix.cc
        RealUpperTriangularMatrix.cc
        RealLowerTriangularMatrix.cc
        ComplexUpperTriangularMatrix.cc
        ComplexLowerTriangularMatrix.cc
        ComplexDiagonalMatrix.cc
        RealBandDiagonalSymmetricMatrix.cc
        ComplexLapackDeterminant.cc
        BandDiagonalHermitianMatrix.cc
        SparseComplexMatrix.cc
        SparseRealMatrix.cc
        LongRationalMatrix.cc
        ComplexUpperHessenbergMatrix.cc
        RealUpperHessenbergMatrix.cc
        IntegerMatrix.cc
        LongIntegerMatrix.cc
)

diagham_add_library(Array
    SOURCES
        Double3DArray.cc
        Complex3DArray.cc
        Complex4DArray.cc
)
