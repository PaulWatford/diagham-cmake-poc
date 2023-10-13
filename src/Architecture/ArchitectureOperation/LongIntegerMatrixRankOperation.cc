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


#include "config.h"
#include "Architecture/ArchitectureOperation/LongIntegerMatrixRankOperation.h"
#include "Matrix/SparseRealMatrix.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include <sys/resource.h>
#include <sys/time.h>


// constructor 
//
// sourceMatrix = pointer to the matrix for which the characteric polynomial should be evaluated
// firstVectorIndex = global index of the first vector that should be processed
// pivotPos = pivot position

LongIntegerMatrixRankOperation::LongIntegerMatrixRankOperation (LongIntegerMatrix* sourceMatrix, int firstVectorIndex, int pivotPos)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceMatrix->GetNbrRow();
  this->SourceMatrix = sourceMatrix;
  this->PivotPos = pivotPos;
  this->FirstVectorIndex = firstVectorIndex;
  this->OperationType = AbstractArchitectureOperation::LongIntegerMatrixRank;
}


// copy constructor 
//
// operation = reference on operation to copy

LongIntegerMatrixRankOperation::LongIntegerMatrixRankOperation(const LongIntegerMatrixRankOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->SourceMatrix = operation.SourceMatrix;
  this->PivotPos = operation.PivotPos;
  this->FirstVectorIndex = operation.FirstVectorIndex;
  this->OperationType = AbstractArchitectureOperation::LongIntegerMatrixRank;
}
  
// destructor
//

LongIntegerMatrixRankOperation::~LongIntegerMatrixRankOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void LongIntegerMatrixRankOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* LongIntegerMatrixRankOperation::Clone()
{
  return new LongIntegerMatrixRankOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool LongIntegerMatrixRankOperation::RawApplyOperation()
{
  int LastComponent = this->FirstComponent + this->NbrComponent;
#ifdef __GMP__
  mpz_t Factor;
  mpz_init(Factor);
#else
  LONGLONG Factor;
#endif
  LongIntegerVector& TmpColumn2 = this->SourceMatrix->Columns[this->FirstVectorIndex - 1];
  for (int TmpFirstNonZero = this->FirstComponent; TmpFirstNonZero < LastComponent; ++TmpFirstNonZero)
    {
#ifdef __GMP__
      mpz_set(Factor, this->SourceMatrix->Columns[TmpFirstNonZero][this->PivotPos]);
#else	      
      Factor = this->SourceMatrix->Columns[TmpFirstNonZero][this->PivotPos];
#endif
      this->SourceMatrix->Columns[TmpFirstNonZero].RescaleAndSubLinearCombination(TmpColumn2[this->PivotPos], Factor, TmpColumn2, this->PivotPos,
										  SourceMatrix->NbrRow - this->PivotPos);
      this->SourceMatrix->Columns[TmpFirstNonZero].Normalize();
    }
  return true;
}

// apply operation for mono processor architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool LongIntegerMatrixRankOperation::ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture)
{
  this->FirstComponent = this->FirstVectorIndex;
  this->NbrComponent = this->SourceMatrix->GetNbrRow() - this->FirstComponent;

  this->RawApplyOperation();

  return true;
}
  
// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool LongIntegerMatrixRankOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  this->FirstComponent = this->FirstVectorIndex;
  this->NbrComponent = this->SourceMatrix->GetNbrRow() - this->FirstComponent;
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  if (Step == 0)
    {
      this->RawApplyOperation();
      return true;
    }
  int TmpFirstComponent = this->FirstVectorIndex;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  LongIntegerMatrixRankOperation** TmpOperations = new LongIntegerMatrixRankOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (LongIntegerMatrixRankOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
    }
   TmpOperations[ReducedNbrThreads] = (LongIntegerMatrixRankOperation*) this->Clone();
   TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->SourceMatrix->GetNbrRow() - TmpFirstComponent);  
   architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);

   architecture->SendJobs();
      

   for (int i = 0; i < architecture->GetNbrThreads(); ++i)
     {
       delete TmpOperations[i];
     }
   delete[] TmpOperations;
   return true;
}


