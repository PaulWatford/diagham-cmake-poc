////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of operation to compute the rank of along integer matrix        //
//                                                                            //
//                        last modification : 12/10/2023                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef LONGINTEGERMATRIXRANKOPERATION_H
#define LONGINTEGERMATRIXRANKOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/LongIntegerMatrix.h"


class LongIntegerMatrixRankOperation: public AbstractArchitectureOperation
{

 protected:

  // pivot position
  int PivotPos;

  // global index of the first vector that should be processed
  int FirstVectorIndex;
  
  // index of the first vector that should be processed locally
  int FirstComponent;
  // number of vectors to process that should be processed locally
  int NbrComponent;

  // pointer to the matrix
  LongIntegerMatrix* SourceMatrix;

  
 public:
  
  // constructor 
  //
  // sourceMatrix = pointer to the matrix for which the characteric polynomial should be evaluated
  // firstVectorIndex = global index of the first vector that should be processed
  // pivotPos = pivot position
  LongIntegerMatrixRankOperation(LongIntegerMatrix* sourceMatrix, int firstVectorIndex, int pivotPos);

  // copy constructor 
  //
  // operation = reference on operation to copy
  LongIntegerMatrixRankOperation(const LongIntegerMatrixRankOperation& operation);
  
  // destructor
  //
  ~LongIntegerMatrixRankOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
    
 protected:

  // apply operation for mono processor architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture);
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // common code for the algorithm initialization
  //
  // trace = pointer to the trace
#ifdef __GMP__
  void AlgorithmInitialization(mpz_t* trace);
#else
  void AlgorithmInitialization(LONGLONG* trace);
#endif
  
};


#endif
