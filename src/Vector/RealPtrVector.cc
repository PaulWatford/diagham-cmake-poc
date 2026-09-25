////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//   class for n dimensional vector each entry pointing to a real number       //
//                                                                            //
//                        last modification : 04/01/2001                      //
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

#include "RealPtrVector.h"
#include "Vector/RealVector.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <fstream>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//

RealPtrVector::RealPtrVector()
{
  this->VectorType = Vector::RealPtrDatas;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->Components = 0;
  this->VectorId = 0;
}

// constructor for an empty real vector (all coordinates set to zero)
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

RealPtrVector::RealPtrVector(int size, bool zeroFlag)
{
  this->VectorType = Vector::RealPtrDatas;
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Components = new double*[this->Dimension + 1]; 
  this->Flag.Initialize();
  this->VectorId = 0;
  if (zeroFlag == true)
    for (int i = 0; i < this->Dimension; i++)
      {
	this->Components[i] = NULL;
      }
}

// copy constructor
//
// vector = vector to copy
// DuplicateFlag = true if datas have to be duplicated

RealPtrVector::RealPtrVector(const RealPtrVector& vector, bool duplicateFlag)
{
  this->VectorType = Vector::RealPtrDatas;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  if (duplicateFlag == false)
    {
      this->Components = vector.Components;
      this->Flag = vector.Flag;
    }
  else
    {
      if (vector.Dimension > 0)
	{
	  this->Flag.Initialize();
	  this->Components = new double*[this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    this->Components[i] = vector.Components[i];
	}
      else
	{
	  this->Components = 0;
	}
    }
}


// copy constructor from a vector (duplicate datas if necessary)
//
// vector = vector to copy

RealPtrVector::RealPtrVector(const Vector& vector)
{
  this->VectorType = Vector::RealPtrDatas;
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->VectorId = vector.VectorId;
  if (vector.VectorType == Vector::RealPtrDatas)
    {
      this->VectorType = Vector::RealPtrDatas;
      this->Components = ((RealPtrVector&) vector).Components;
      this->Flag = ((RealPtrVector&) vector).Flag;
    }
  else
      {
	this->Components = 0;
	this->Flag.Initialize();
      }
}

#ifdef __MPI__

// constructor from informations sent using MPI ... not checked, probably buggy!
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the vector
// broadcast = true if the vector is broadcasted

RealPtrVector::RealPtrVector(const MPI_Comm& communicator, int id, bool broadcast)
{
  this->VectorType = Vector::RealPtrDatas;
  int TmpArray[3];
  if (broadcast == true)
    {
      MPI_Bcast(TmpArray, 3, MPI_INT, id, communicator);
    }
  else
    {
      MPI_Status TmpMPIStatus;
      MPI_Recv(TmpArray, 3, MPI_INT, id, 1, communicator, &TmpMPIStatus);
    }
  this->Dimension = TmpArray[0];
  this->VectorId = TmpArray[1];
  this->Components = new double *[this->Dimension + 1];
  if (TmpArray[2] == 1)
    {
      for (int i = 0; i <= this->Dimension; ++i) 
	this->Components[i] = NULL;
    }
  else
    {
      if (TmpArray[2] == 2)
	{
	  if (broadcast == true)
	    {
	      MPI_Bcast(this->Components, this->Dimension, MPI_DOUBLE, id, communicator);
	    }
	  else
	    {
	      MPI_Status TmpMPIStatus;
	      MPI_Recv(this->Components, this->Dimension, MPI_DOUBLE, id, 1, communicator, &TmpMPIStatus);
	    }
	}
    }
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
}

#endif

// destructor
//

RealPtrVector::~RealPtrVector ()
{
  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

RealPtrVector& RealPtrVector::operator = (const RealPtrVector& vector)
{
  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Flag = vector.Flag;
  this->VectorId = vector.VectorId;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.Dimension;
  return *this;
}


// assignement from a vector (duplicate datas if necessary)
//
// vector = vector to assign
// return value = reference on current vector

RealPtrVector& RealPtrVector::operator = (const Vector& vector)
{
  if (vector.VectorType == Vector::RealPtrDatas)
    {
      return ((*this) = (RealPtrVector&) vector);
    }
  return *this;
}

// Resize vector
//
// dimension = new dimension

void RealPtrVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  double** TmpVector = new double*[dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    TmpVector[i] = this->Components[i];
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->Components = TmpVector;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void RealPtrVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  double** TmpVector = new double *[dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    TmpVector[i] = this->Components[i];
  for (int i = this->Dimension; i < dimension; i++)
    TmpVector[i] = NULL;  
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->Components = TmpVector;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

RealPtrVector& RealPtrVector::Copy (RealPtrVector& vector)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  this->Localize();
  vector.Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] = vector.Components[i];
  this->Delocalize();
  vector.Delocalize();
  return *this;
}


// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

double operator * (RealPtrVector& V1, RealPtrVector& V2)
{
  V1.Localize();
  V2.Localize();
  double **MyComponentPtr=V1.Components;
  double **OtherComponentPtr=V2.Components;
  double x = **MyComponentPtr++ * **OtherComponentPtr++;
  int L=V1.Dimension;
  for (int i = 1; i < L; ++i)
     x += **MyComponentPtr++ * **OtherComponentPtr++;  
  V1.Delocalize();
  V2.Delocalize();
  return x;
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

double operator * (RealPtrVector& V1, RealVector& V2)
{
  V1.Localize();
  V2.Localize();
  double **MyComponentPtr=V1.Components;
  double *OtherComponentPtr=V2.Components;
  double x = **MyComponentPtr++ * *OtherComponentPtr++;
  int L=V1.Dimension;
  //#pragma simd
  for (int i = 1; i < L; ++i)
     x += **MyComponentPtr++ * *OtherComponentPtr++;     
  V1.Delocalize();
  V2.Delocalize();
  return x;
}


// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

double operator * (RealVector& V1, RealPtrVector& V2)
{
  V1.Localize();
  V2.Localize();  
  double *MyComponentPtr=V1.Components;
  double **OtherComponentPtr=V2.Components;
  int L=V1.Dimension;
  //#define OTHER_METHOD
#ifndef OTHER_METHOD
  // baseline method: use pointers
  double x = *MyComponentPtr++ * **OtherComponentPtr++;
  // this pragma works, but does not appear to affect computer time significantly
  //# pragma simd
  for (int i = 1; i < L; ++i)
  // double *end = V1.Components+V1.Dimension;
  // for (; MyComponentPtr < end; )
    x += *MyComponentPtr++ * **OtherComponentPtr++;
#else
  // alternative method - no pointer arithmetic, but dereferencing: 2x slower than above!
  double x=0.0;
  // this pragma does not lead to vectorization:
#pragma simd
  for (int i = 0; i < L; ++i)
    {
      x += V1[i] * *V2[i];
    }
#endif
  V1.Delocalize();
  V2.Delocalize();
  return x;
}

  
// Output Stream overload
//
// str = reference on output stream
// v = vector to print
// return value = reference on output stream

ostream& operator << (ostream& str, RealPtrVector& v)
{
  v.Localize();
  for (int i = 0; i < v.Dimension; ++i)
    {
      str << *(v.Components[i]) << endl;
    }
  v.Delocalize();
  return str;
}


#ifdef __MPI__

// send a vector to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& RealPtrVector::SendVector(const MPI_Comm& communicator, int id)
{
//   MPI_Send(&this->VectorType, 1, MPI_INT, id, 1);
//   MPI_Send(&this->Dimension, 1, MPI_INT, id, 1); 
//   int Acknowledge = 0;
//   MPI_Recv(&Acknowledge, 1, MPI_INT, id, 1);
//   if (Acknowledge != 0)
//     return *this;
//   MPI_Send(this->Components, this->Dimension, MPI_DOUBLE, id, 1); 
  return *this;
}

// broadcast a vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// return value = reference on the current vector

Vector& RealPtrVector::BroadcastVector(const MPI_Comm& communicator,  int id)
{
//   int TmpVectorType = this->VectorType;
//   int TmpDimension = this->Dimension;
//   int Acknowledge = 0;
//   MPI_Bcast(&TmpVectorType, 1, MPI_INT, id);
//   MPI_Bcast(&TmpDimension, 1, MPI_INT, id);
//   if (this->VectorType != TmpVectorType)
//     {
//       Acknowledge = 1;
//     }
//   if (id != MPI_Get_rank())
//     MPI_Send(&Acknowledge, 1, MPI_INT, id, 1);      
//   else
//     {
//       int NbrMPINodes = MPI_Get_size();
//       bool Flag = false;
//       for (int i = 0; i < NbrMPINodes; ++i)
// 	if (id != i)
// 	  {
// 	    MPI_Recv(&Acknowledge, 1, MPI_INT, i, 1);      
// 	    if (Acknowledge == 1)
// 	      Flag = true;
// 	  }
//       if (Flag == true)
// 	Acknowledge = 1;
//     }
//   MPI_Bcast(&Acknowledge, 1, MPI_INT, id);
//   if (Acknowledge != 0)
//     return *this;
//   if (TmpDimension != this->Dimension)
//     {
//       this->Resize(TmpDimension);      
//     }
//   MPI_Bcast(this->Components, this->Dimension, MPI_DOUBLE, id);
  return *this;
}

// broadcast part of vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
// nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
// return value = reference on the current vector

Vector& RealPtrVector::BroadcastPartialVector(const MPI_Comm& communicator, int id, int firstComponent, int nbrComponent)
{
//   int TmpVectorType = this->VectorType;
//   int TmpDimension = this->Dimension;
//   int Acknowledge = 0;
//   MPI_Bcast(&TmpVectorType, 1, MPI_INT, id);
//   MPI_Bcast(&TmpDimension, 1, MPI_INT, id);
//   MPI_Bcast(&firstComponent, 1, MPI_INT, id);
//   MPI_Bcast(&nbrComponent, 1, MPI_INT, id);
//   if (this->VectorType != TmpVectorType)
//     {
//       Acknowledge = 1;
//     }
//   if (id != MPI_Get_rank())
//     MPI_Send(&Acknowledge, 1, MPI_INT, id, 1);      
//   else
//     {
//       int NbrMPINodes = MPI_Get_size();
//       bool Flag = false;
//       for (int i = 0; i < NbrMPINodes; ++i)
// 	if (id != i)
// 	  {
// 	    MPI_Recv(&Acknowledge, 1, MPI_INT, i, 1);      
// 	    if (Acknowledge == 1)
// 	      Flag = true;
// 	  }
//       if (Flag == true)
// 	Acknowledge = 1;
//     }
//   MPI_Bcast(&Acknowledge, 1, MPI_INT, id);
//   if (Acknowledge != 0)
//     return *this;
//   if (TmpDimension != this->Dimension)
//     {
//       this->Resize(TmpDimension);      
//     }
//   MPI_Bcast(this->Components + firstComponent, nbrComponent, MPI_DOUBLE, id);
  return *this;
}

// receive a vector from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current vector

Vector& RealPtrVector::ReceiveVector(const MPI_Comm& communicator, int id)
{
  int TmpVectorType = 0;
//   int TmpDimension = 0;
//   MPI_Recv(&TmpVectorType, 1, MPI_INT, id, 1);
//   MPI_Recv(&TmpDimension, 1, MPI_INT, id, 1); 
//   if (TmpVectorType != this->VectorType)
//     {
//       TmpDimension = 1;
//       MPI_Send(&TmpDimension, 1, MPI_INT, id, 1);
//       return *this;
//     }
//   else
//     {
//       if (TmpDimension != this->Dimension)
// 	{
// 	  this->Resize(TmpDimension);      
// 	}
//       TmpDimension = 0;
//       MPI_Send(&TmpDimension, 1, MPI_INT, id, 1);
//     }
//   MPI_Recv(this->Components, this->Dimension, MPI_DOUBLE, id, 1); 
  return *this;
}

// add current vector to the current vector of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& RealPtrVector::SumVector(const MPI_Comm& communicator, int id)
{
//   int TmpVectorType = this->VectorType;
//   int TmpDimension = this->Dimension;
//   int Acknowledge = 0;
//   MPI_Bcast(&TmpVectorType, 1, MPI_INT, id);
//   MPI_Bcast(&TmpDimension, 1, MPI_INT, id);
//   if ((this->VectorType != TmpVectorType) || (this->Dimension != TmpDimension))
//     {
//       Acknowledge = 1;
//     }
//   if (id != MPI_Get_rank())
//     MPI_Send(&Acknowledge, 1, MPI_INT, id, 1);      
//   else
//     {
//       int NbrMPINodes = MPI_Get_size();
//       bool Flag = false;
//       for (int i = 0; i < NbrMPINodes; ++i)
// 	if (id != i)
// 	  {
// 	    MPI_Recv(&Acknowledge, 1, MPI_INT, i, 1);      
// 	    if (Acknowledge == 1)
// 	      Flag = true;
// 	  }
//       if (Flag == true)
// 	Acknowledge = 1;
//     }
//   MPI_Bcast(&Acknowledge, 1, MPI_INT, id);
//   if (Acknowledge != 0)
//     {
//       return *this;
//     }
//   double* TmpComponents = 0;
//   if (id == MPI_Get_rank())
//     {
//       TmpComponents = new double [this->Dimension];
//     }
//   MPI_Reduce(this->Components, TmpComponents, this->Dimension, MPI_DOUBLE, MPI_SUM, id);
//   if (id == MPI_Get_rank())
//     {
//       for (int i = 0; i < this->Dimension; ++i)
// 	this->Components[i] = TmpComponents[i];
//       delete[] TmpComponents;
//     }
  return *this;
}

// create a new vector on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* RealPtrVector::BroadcastClone(const MPI_Comm& communicator, int id)
{
//   if (id == MPI_Get_rank())
//     {
//       MPI_Bcast(&this->VectorType, 1, MPI_INT, id);
//       int TmpArray[3];
//       TmpArray[0] = this->Dimension;
//       TmpArray[1] = this->VectorId;
//       TmpArray[2] = 2;
//       MPI_Bcast(TmpArray, 3, MPI_INT, id);      
//       MPI_Bcast(this->Components, this->Dimension, MPI_DOUBLE, id);      
//     }
//   else
//     {
//       int Type = 0;
//       MPI_Bcast(&Type, 1, MPI_INT, id);  
//       return new RealPtrVector(communicator, id);
//     }
  return 0;
}

// create a new vector on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* RealPtrVector::BroadcastEmptyClone(const MPI_Comm& communicator, int id, bool zeroFlag)
{
//   if (id == MPI_Get_rank())
//     {
//       MPI_Bcast(&this->VectorType, 1, MPI_INT, id);
//       int TmpArray[3];
//       TmpArray[0] = this->Dimension;
//       TmpArray[1] = this->VectorId;
//       TmpArray[2] = 0;
//       if (zeroFlag == true)
// 	{
// 	  TmpArray[2] = 1;
// 	}
//       MPI_Bcast(TmpArray, 3, MPI_INT, id);      
//     }
//   else
//     {
//       int Type = 0;
//       MPI_Bcast(&Type, 1, MPI_INT, id);  
//       return new RealPtrVector(communicator, id);
//     }
  return 0;
}

#endif
