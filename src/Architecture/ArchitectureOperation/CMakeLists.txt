# Auto-generated CMakeLists.txt for src/Architecture/ArchitectureOperation
# Regenerate with: python3 scripts_cmake/extract_autotools.py

diagham_add_library(ArchitectureOperation
    SOURCES
        AbstractArchitectureOperation.cc
        VectorHamiltonianMultiplyOperation.cc
        VectorOperatorMultiplyOperation.cc
        MultipleVectorHamiltonianMultiplyOperation.cc
        MultipleRealScalarProductOperation.cc
        AddRealLinearCombinationOperation.cc
        MultipleComplexScalarProductOperation.cc
        AddComplexLinearCombinationOperation.cc
        MatrixMatrixMultiplyOperation.cc
        AbstractPrecalculationOperation.cc
        AbstractScalarSumOperation.cc
        MainTaskOperation.cc
        AbstractArchitectureOperationManager.cc
        ArchitectureBaseOperationManager.cc
        OperatorMatrixElementOperation.cc
        HamiltonianFullDiagonalizeOperation.cc
        SparseMatrixMatrixMultiplyOperation.cc
        GenericHamiltonianPrecalculationOperation.cc
        MultipleOperatorMatrixElementOperation.cc
        VectorTensorMultiplicationCoreOperation.cc
        VectorSparseTensorMultiplyOperation.cc
        TensorMatrixContractionOperation.cc
        MPOApplyOnTensorOperation.cc
        TensorVectorContractionOperation.cc
        OperatorMultipleMatrixElementOperation.cc
        MatrixFullDiagonalizeOperation.cc
        MultipleVectorOperatorMultiplyOperation.cc
        HermitianMatrixFromMatrixOperation.cc
        LongIntegerMatrixCharacteristicPolynomialOperation.cc
        LongIntegerMatrixRankOperation.cc
)
