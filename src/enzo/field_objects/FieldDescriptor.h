//
// Field Descriptor Object
//
// Authors: Matthew Turk
//          Greg Bryan

#ifndef __FIELD_DESCRIPTOR_H__
#define __FIELD_DESCRIPTOR_H__

class Grid;

class FieldDescriptor 
{

    public:
      // Constructors
      FieldDescriptor();
      FieldDescriptor(FieldDescriptor* BaseDefinition,
                      int CellDimensions[MAX_DIMENSIONS],
                      long_int LeftEdge[MAX_DIMENSIONS],
                      float **FieldPointer = NULL,
                      int SkipValueAllocation = 0);
      FieldDescriptor(CenteringType ValueCentering, int Rank,
                      int CellDimensions[MAX_DIMENSIONS],
                      long_int LeftEdge[MAX_DIMENSIONS],
                      InterpolationType InterpolationMethod,
                      const char* Name,
                      const char* UnitsName,
                      float **FieldPointer = NULL);

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
      float *GetValues();
      int GetSize();
      void GetCellDimensions(int Dimensions[MAX_DIMENSIONS]);
      void GetFieldDimensions(int Dimensions[MAX_DIMENSIONS]);
      void GetLeftEdge(long_int LeftEdge[MAX_DIMENSIONS]);
      InterpolationType GetInterpolationMethod();

      void GetOverlapRegion(FieldDescriptor *Other,
                int LeftEdgeThis[MAX_DIMENSIONS],
                int LeftEdgeOther[MAX_DIMENSIONS],
                int CopyDims[MAX_DIMENSIONS]);

      // Mathematical Operations

      float Min();
      float Min(int *LeftEdge, int *RightEdge);
      float Max();
      float Max(int *LeftEdge, int *RightEdge);
      float Sum();
      float Sum(int *LeftEdge, int *RightEdge);

      // Operations from other FieldDescriptors

      void CopyFrom(FieldDescriptor *Other);
      void CopyFrom(float val);

      void Add(FieldDescriptor *Other);
      void Add(float val);

      void Subtract(FieldDescriptor *Other);
      void Subtract(float val);

      void Multiply(FieldDescriptor *Other);
      void Multiply(float val);

      void Divide(FieldDescriptor *Other);
      void Divide(float val);

      // Now the fun stuff: operator overloading!
      // We'll do in-place operators first.

      FieldDescriptor &operator+=(FieldDescriptor *Other);
      FieldDescriptor &operator-=(FieldDescriptor *Other);
      FieldDescriptor &operator*=(FieldDescriptor *Other);
      FieldDescriptor &operator/=(FieldDescriptor *Other);

      FieldDescriptor &operator+=(float val);
      FieldDescriptor &operator-=(float val);
      FieldDescriptor &operator*=(float val);
      FieldDescriptor &operator/=(float val);

      int CanCombine(FieldDescriptor *Other);

    protected: 

      template <MathFunction function>
        float UnaryAccumulator(
            int LeftEdge[MAX_DIMENSIONS],
            int RightEdge[MAX_DIMENSIONS],
            float InitialValue);

      template <MathFunction function>
      void InPlaceBinaryOperation(
          FieldDescriptor *Other,
          int LeftEdgeThis[MAX_DIMENSIONS],
          int LeftEdgeOther[MAX_DIMENSIONS],
          int CopyDims[MAX_DIMENSIONS]);

      template <MathFunction function>
      void InPlaceBinaryOperation(
          float OtherValue,
          int LeftEdge[MAX_DIMENSIONS],
          int CopyDims[MAX_DIMENSIONS]);

      void DeallocateIfNeeded();
      void SetPointer(float **NewPointer, int SkipValueAllocation = 0);
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
      long_int LeftEdge[MAX_DIMENSIONS];
      float **FieldPointer;
      int DeallocateFieldPointer;
      int DeallocateFieldValues;

};

#endif
