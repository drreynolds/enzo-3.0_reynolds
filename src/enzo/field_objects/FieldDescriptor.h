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
      //
      // This is the constructor for a bare field descriptor.  It is unlikely
      // to be used in practice, as it will then require setting up the
      // dimensions and so on at a later time, which is non-trivial.  However,
      // it is often used internally for copying and so on.
      FieldDescriptor();
      // This constructor starts from an *existing* field definition (see
      // FieldRegistry.C for existing defined fields) and then constructs that
      // field (and its corresponding centering, etc) in a new object, which
      // can then be modified.  This is the most common constructor.
      //
      // This function accepts:
      //    CellDimensions    
      //        This should be the number of cell dimensions, not the number of
      //        values expected for the specified centering type.  Internally,
      //        the field object will calculate the necessary number of values
      //        for the specified centering.
      //    LeftEdge
      //        This is the left edge, in (arbitrary, positive-definite)
      //        integers of the region of space that this field occupies.
      //        Often this is calculated by taking the left edge of a grid
      //        object and dividing by the local cell width; in this way, it's
      //        an integer coordinate constructed at the local refinement
      //        level.
      //    FieldPointer
      //        Either a pointer to a field, or a pointer to a pointer to a
      //        field.  Can be NULL, or a pointer to a NULL field.  When used
      //        with Enzo grid objects, typically a FieldDescriptor will be
      //        instantiated with a pointer to a BaryonField pointer.  In this
      //        way, it will be affiliated with a NULL pointer, or a pointer to
      //        a real field, and it will be able to live outside of any
      //        allocations and deallocations that occur.  If a pointer is
      //        supplied, it will be "owned" by this object.
      //    SkipValueAllocation
      //        This specifies whether we should skip new-ing a field if the
      //        field pointer is NULL.  Typically, when working with Enzo grid
      //        objects, we will skip value allocation.
      FieldDescriptor(FieldDescriptor* BaseDefinition,
                      int CellDimensions[MAX_DIMENSIONS],
                      long_int LeftEdge[MAX_DIMENSIONS],
                      float **FieldPointer = NULL,
                      int SkipValueAllocation = 0);
      // This constructor is the most raw method of constructing a
      // FieldDescriptor, and it allows fields to be created from no previous
      // information.  It may be useful for one-off fields that do not persist
      // for a long time, such as simple derived values and the like.
      //
      // This function accepts:
      //    CenteringType
      //        This describes the type of centering the field values have.
      //        The full list can be found in FieldDefinitions.h, but it
      //        includes things like face, edge, corner, cell.
      //    Rank
      //        The rank of the field, i.e., the dimensionality.  1, 2, 3
      //        typically.
      //    CellDimensions    
      //        This should be the number of cell dimensions, not the number of
      //        values expected for the specified centering type.  Internally,
      //        the field object will calculate the necessary number of values
      //        for the specified centering.
      //    LeftEdge
      //        This is the left edge, in (arbitrary, positive-definite)
      //        integers of the region of space that this field occupies.
      //        Often this is calculated by taking the left edge of a grid
      //        object and dividing by the local cell width; in this way, it's
      //        an integer coordinate constructed at the local refinement
      //        level.
      //    InterpolationMethod
      //        When we are interpolating between levels, do we need to do
      //        anything special during the process?  Options include
      //        multiplying by density, not interpolating at all, directly
      //        interpolating, or undefined.  New methods can be easily added.
      //    Name
      //        The name of the field.  Will be used for fast lookups.
      //    UnitsName
      //        The name of the units for this field.  Not currently used.
      //    FieldPointer
      //        Either a pointer to a field, or a pointer to a pointer to a
      //        field.  Can be NULL, or a pointer to a NULL field.  When used
      //        with Enzo grid objects, typically a FieldDescriptor will be
      //        instantiated with a pointer to a BaryonField pointer.  In this
      //        way, it will be affiliated with a NULL pointer, or a pointer to
      //        a real field, and it will be able to live outside of any
      //        allocations and deallocations that occur.  If a pointer is
      //        supplied, it will be "owned" by this object.  If NULL, it will
      //        be allocated.
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

      // This function, typically only used internally but still public, is
      // designed to determine the overlap in *cell* space of two different
      // FieldDescriptors.  Note that any other functions that allow left edge
      // and copydims to be specified will internally calculate the offset and
      // necessary indices to account for centering of values.
      //
      //    FieldDescriptor
      //        The other field descriptor, with which the calculated overlap
      //        is desired.
      //    LeftEdgeThis
      //        This (output parameter) is the left edge in *local* coordinates
      //        of *cells* for this field descriptor for the overlap between
      //        the two field descriptors.
      //    LeftEdgeOther
      //        This (output parameter) is the left edge in *local* coordinates
      //        of *cells* for the *other* field descriptor.for the overlap
      //        between the two field descriptors.
      //    CopyDims
      //        The number of *cell* values that overlap.
      void GetOverlapRegion(FieldDescriptor *Other,
                int LeftEdgeThis[MAX_DIMENSIONS],
                int LeftEdgeOther[MAX_DIMENSIONS],
                int CopyDims[MAX_DIMENSIONS]);

      // Mathematical Operations

      // These operations return a single value.

      float Min();
      float Min(int *LeftEdge, int *RightEdge);
      float Max();
      float Max(int *LeftEdge, int *RightEdge);
      float Sum();
      float Sum(int *LeftEdge, int *RightEdge);

      // Operations from other FieldDescriptors

      // These operations accept either another FieldDescriptor (in which case
      // the correct bounds for the computation will be computed) or a scalar
      // value.

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
      
      // These are in-place operators that modify the values of the current
      // field.

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

      // These "UnaryAccumulator" functions are templated functions that enable
      // a function to collects a single value to be returned at the end of
      // computation.  These could, for instance, be things like summations,
      // standard deviations, etc.  Note that the left edge and right edge
      // values are integers in local *cell* coordinates, and computations that
      // are conducted on face and edge centered values will compute internally
      // the correct indices.
      template <MathFunction function>
        float UnaryAccumulator(
            int LeftEdge[MAX_DIMENSIONS],
            int RightEdge[MAX_DIMENSIONS],
            float InitialValue);

      // These "InPlaceBinary" functions are templated functions that enable
      // a function to process two different field descriptors and modify
      // *this* one in place.  This could, for instance, be multiplying two
      // field descriptors or subtracting them.  As above, the dimensions are
      // specified in cell coordinates and will be converted to indices in the
      // fields themselves.
      template <MathFunction function>
      void InPlaceBinaryOperation(
          FieldDescriptor *Other,
          int LeftEdgeThis[MAX_DIMENSIONS],
          int LeftEdgeOther[MAX_DIMENSIONS],
          int CopyDims[MAX_DIMENSIONS]);

      // This inplace binary operation enables a single scalar value to be used
      // as the other value for a binary operation.  For instance, this could
      // be multiplying by a constant.
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
