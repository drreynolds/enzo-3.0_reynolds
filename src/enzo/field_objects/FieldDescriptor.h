//
// Field Descriptor Object
//
// Authors: Matthew Turk
//          Greg Bryan

#ifndef __FIELD_DESCRIPTOR_H__
#define __FIELD_DESCRIPTOR_H__

#include <string.h>
#include <string>

class Grid;

class FieldDescriptor 
{

    public:
      // Constructors
      FieldDescriptor();
      FieldDescriptor(FieldDescriptor* BaseDefinition,
                      int CellDimensions[MAX_DIMENSIONS],
                      long long LeftEdge[MAX_DIMENSIONS],
                      double **FieldPointer = NULL);
      FieldDescriptor(CenteringType ValueCentering, int Rank,
                      int CellDimensions[MAX_DIMENSIONS],
                      long long LeftEdge[MAX_DIMENSIONS],
                      InterpolationType InterpolationMethod,
                      const char* Name,
                      const char* UnitsName,
                      double **FieldPointer = NULL);

      // Destructor
      ~FieldDescriptor();

      // Setters and Getters

      FieldDescriptor *Duplicate(char *NewName = NULL, char *NewUnitsName = NULL);

      int Index(int i, int j, int k);
      CenteringType GetValueCentering();
      char* GetName();
      void SetName(std::string NewName);
      void SetName(const char* NewName);
      char* GetUnitsName();
      double *GetValues();
      int GetSize();
      void GetCellDimensions(int Dimensions[MAX_DIMENSIONS]);
      void GetFieldDimensions(int Dimensions[MAX_DIMENSIONS]);
      void GetLeftEdge(long long LeftEdge[MAX_DIMENSIONS]);

      void GetOverlapRegion(FieldDescriptor *Other,
                int LeftEdgeThis[MAX_DIMENSIONS],
                int LeftEdgeOther[MAX_DIMENSIONS],
                int CopyDims[MAX_DIMENSIONS]);

      // Mathematical Operations

      double min();
      double min(int *LeftEdge, int *RightEdge);
      double max();
      double max(int *LeftEdge, int *RightEdge);
      double sum();
      double sum(int *LeftEdge, int *RightEdge);

      // Operations from other FieldDescriptors

      void CopyFrom(FieldDescriptor *Other);
      void CopyFrom(double val);

      void Add(FieldDescriptor *Other);
      void Add(double val);

      void Subtract(FieldDescriptor *Other);
      void Subtract(double val);

      void Multiply(FieldDescriptor *Other);
      void Multiply(double val);

      void Divide(FieldDescriptor *Other);
      void Divide(double val);

      // Now the fun stuff: operator overloading!
      // We'll do in-place operators first.

      FieldDescriptor &operator+=(FieldDescriptor *Other);
      FieldDescriptor &operator-=(FieldDescriptor *Other);
      FieldDescriptor &operator*=(FieldDescriptor *Other);
      FieldDescriptor &operator/=(FieldDescriptor *Other);

      FieldDescriptor &operator+=(double val);
      FieldDescriptor &operator-=(double val);
      FieldDescriptor &operator*=(double val);
      FieldDescriptor &operator/=(double val);

      int CanCombine(FieldDescriptor *Other);

    protected: 

      template <MathFunction function>
        double UnaryAccumulator(
            int LeftEdge[MAX_DIMENSIONS],
            int RightEdge[MAX_DIMENSIONS],
            double InitialValue);

      template <MathFunction function>
      void InPlaceBinaryOperation(
          FieldDescriptor *Other,
          int LeftEdgeThis[MAX_DIMENSIONS],
          int LeftEdgeOther[MAX_DIMENSIONS],
          int CopyDims[MAX_DIMENSIONS]);

      template <MathFunction function>
      void InPlaceBinaryOperation(
          double OtherValue,
          int LeftEdge[MAX_DIMENSIONS],
          int CopyDims[MAX_DIMENSIONS]);

      void DeallocateIfNeeded();
      void SetPointer(double **NewPointer);
      void AllocateFieldValues();
      void AllocateFieldPointer();

      CenteringType ValueCentering;
      int Rank;
      InterpolationType InterpolationMethod;
      void GetFieldExtension(int Extension[MAX_DIMENSIONS]);

    private:
      char* Name;
      char* UnitsName;
      int CellDimensions[MAX_DIMENSIONS];
      long long LeftEdge[MAX_DIMENSIONS];
      double **FieldPointer;
      int DeallocateFieldPointer;
      int DeallocateFieldValues;

};

#endif
