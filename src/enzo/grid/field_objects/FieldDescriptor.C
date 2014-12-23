//
// Field Descriptor Object Implementation
//
// Authors: Matthew Turk
//          Greg Bryan

#include "FieldObjects.h"

FieldDescriptor::FieldDescriptor(
    FieldDescriptor *BaseDefinition,
    int CellDimensions[MAX_DIMENSIONS],
    long_int LeftEdge[MAX_DIMENSIONS],
    float **FieldPointer, int SkipValueAllocation, int Level) {
  int dim;
  this->ValueCentering = BaseDefinition->ValueCentering;
  this->Rank = BaseDefinition->Rank;
  for (dim = 0; dim < this->Rank; dim++) {
    this->CellDimensions[dim] = CellDimensions[dim];
    this->LeftEdge[dim] = LeftEdge[dim];
  }
  this->Level = Level;
  this->InterpolationMethod = BaseDefinition->InterpolationMethod;
  // Do we want to strdup or just assign the pointer?
  this->Name = strdup(BaseDefinition->GetName());
  this->UnitsName = strdup(BaseDefinition->GetUnitsName());

  // This sets up how we handle the values we point to.
  this->DeallocateFieldValues = this->DeallocateFieldPointer = 0;
  this->SetPointer(FieldPointer, SkipValueAllocation);

#ifdef FIELD_DEBUG
  fprintf(stderr, "OWNERSHIP OF %s IS %"ISYM" %"ISYM"\n",
          this->GetName(), this->DeallocateFieldValues,
          this->DeallocateFieldPointer);
#endif
}

FieldDescriptor::FieldDescriptor(
    CenteringType ValueCentering, int Rank,
    int CellDimensions[MAX_DIMENSIONS],
    long_int LeftEdge[MAX_DIMENSIONS],
    InterpolationType InterpolationMethod,
    const char* Name,
    const char* UnitsName,
    float **FieldPointer,
    int Level) {
  int dim;
  this->ValueCentering = ValueCentering;
  this->Rank = Rank;
  for (dim = 0; dim < this->Rank; dim++) {
    this->CellDimensions[dim] = CellDimensions[dim];
    this->LeftEdge[dim] = LeftEdge[dim];
  }
  this->Level = Level;
  this->InterpolationMethod = InterpolationMethod;
  // Do we want to strdup or just assign the pointer?
  this->Name = strdup(Name);
  this->UnitsName = strdup(UnitsName);

  // This sets up how we handle the values we point to.
  this->DeallocateFieldValues = this->DeallocateFieldPointer = 0;
  this->SetPointer(FieldPointer);
}

FieldDescriptor::~FieldDescriptor() {
  // If we use strdup then we need to de-allocate
  this->DeallocateIfNeeded();
  free(this->Name);
  free(this->UnitsName);
};

void FieldDescriptor::DeallocateIfNeeded() {
  if ((this->DeallocateFieldValues == 1) &&
      (this->FieldPointer[0] != NULL)) {
    delete this->FieldPointer[0];
    this->FieldPointer[0] = NULL;
  }
  if ((this->DeallocateFieldPointer == 1) &&
      (this->FieldPointer != NULL)) {
    delete this->FieldPointer;
    this->FieldPointer = NULL;
  }
}

void FieldDescriptor::AllocateFieldValues() {
  // We do *not* want this to get lost on stack deallocation
  float *v = new float[this->GetSize()];
  this->FieldPointer[0] = v;
  this->DeallocateFieldValues = 1;
}

void FieldDescriptor::AllocateFieldPointer() {
  this->FieldPointer = new float*[1];
  this->FieldPointer[0] = NULL;
  this->DeallocateFieldPointer = 1;
}

FieldDescriptor* FieldDescriptor::Duplicate(char *NewName, char *NewUnitsName) {
  int i;
  float *vn;
  if (NewName == NULL) NewName = this->Name;
  if (NewUnitsName == NULL) NewUnitsName = this->UnitsName;
  // Note that this will allocate and set the bits to deallocate later
  FieldDescriptor *fd = new FieldDescriptor(this->ValueCentering, this->Rank,
        this->CellDimensions, this->LeftEdge, this->InterpolationMethod,
        NewName, NewUnitsName, NULL, this->Level);
  if(this->FieldPointer != NULL) {
    // We copy if there are any FieldPointer
    float *vo = fd->GetValues();
    float *vt = this->GetValues();
    for (i = 0; i < this->GetSize(); i++) {
        vo[i] = vt[i];
    }
  }
  return fd;
}

int FieldDescriptor::Index(int i, int j, int k) {
  /* This should be inlined or sped up or something. */
  int FieldDimensions[MAX_DIMENSIONS];
  this->GetFieldDimensions(FieldDimensions);
  return ((k*FieldDimensions[1]) + j)*FieldDimensions[0] + i;
}

CenteringType FieldDescriptor::GetValueCentering() {
  return this->ValueCentering;
}

char* FieldDescriptor::GetName() {
  return this->Name;
};

char* FieldDescriptor::GetUnitsName() {
  return this->Name;
};

void FieldDescriptor::SetName(std::string NewName) {
  // We own the memory
  if (this->Name != NULL) free(this->Name);
  this->Name = strdup(NewName.c_str());
}

void FieldDescriptor::SetName(const char* NewName) {
  // We own the memory
  if (this->Name != NULL) free(this->Name);
  this->Name = strdup(NewName);
}

float *FieldDescriptor::GetValues() {
  if(this->FieldPointer == NULL) return NULL;
  return *(this->FieldPointer);
}

void FieldDescriptor::SetPointer(float **NewPointer, int SkipValueAllocation) {
  if (NewPointer == NULL) {
    this->AllocateFieldPointer();
    NewPointer = this->FieldPointer;
  } else {
    this->FieldPointer = NewPointer;
  }
  // If we need to, also allocate values
  if ((NewPointer[0] == NULL) && (SkipValueAllocation == 0)) {
    this->AllocateFieldValues();
  }
}

int FieldDescriptor::GetSize() {
  int dim, Dimensions[MAX_DIMENSIONS], size = 1;
  this->GetFieldDimensions(Dimensions);
  for (dim = 0; dim < this->Rank; dim++) {
    size *= Dimensions[dim];
  }
  return size;
}

void FieldDescriptor::GetCellDimensions(int Dimensions[MAX_DIMENSIONS]) {
  int dim;
  for (dim = 0; dim < this->Rank; dim++) {
    Dimensions[dim] = this->CellDimensions[dim];
  }
}

void FieldDescriptor::GetFieldDimensions(int Dimensions[MAX_DIMENSIONS]) {
  int dim;
  int Extension[MAX_DIMENSIONS];
  this->GetCellDimensions(Dimensions);
  this->GetFieldExtension(Extension);
  for (dim = 0; dim < this->Rank; dim++) {
    Dimensions[dim] += Extension[dim];
  }
}

void FieldDescriptor::GetFieldExtension(int Extension[MAX_DIMENSIONS]) {
  switch (this->ValueCentering) {
    case CellCentered:
      Extension[0] = 0;
      Extension[1] = 0;
      Extension[2] = 0;
      break;
    case CornerCentered:
      Extension[0] = 1;
      Extension[1] = 1;
      Extension[2] = 1;
      break;
    case FaceCenteredX:
      Extension[0] = 1;
      Extension[1] = 0;
      Extension[2] = 0;
      break;
    case FaceCenteredY:
      Extension[0] = 0;
      Extension[1] = 1;
      Extension[2] = 0;
      break;
    case FaceCenteredZ:
      Extension[0] = 0;
      Extension[1] = 0;
      Extension[2] = 1;
      break;
    case EdgeCenteredX:
      Extension[0] = 0;
      Extension[1] = 1;
      Extension[2] = 1;
      break;
    case EdgeCenteredY:
      Extension[0] = 1;
      Extension[1] = 0;
      Extension[2] = 1;
      break;
    case EdgeCenteredZ:
      Extension[0] = 1;
      Extension[1] = 1;
      Extension[2] = 0;
      break;
  }
}

void FieldDescriptor::GetLeftEdge(long_int LeftEdge[MAX_DIMENSIONS]) {
    int dim;
    for (dim = 0; dim < this->Rank; dim++) {
        LeftEdge[dim] = this->LeftEdge[dim];
    }
}

InterpolationType FieldDescriptor::GetInterpolationMethod(){
    return this->InterpolationMethod;
}

void FieldDescriptor::GetOverlapRegion(FieldDescriptor *Other,
    int LeftEdgeThis[MAX_DIMENSIONS], // These are out values
    int LeftEdgeOther[MAX_DIMENSIONS],
    int CopyDims[MAX_DIMENSIONS]) {
      // Only CelL overlap
      int OtherDims[MAX_DIMENSIONS];
      long_int OverlapStart, OverlapEnd;
      long_int GlobalLeftEdgeThis[MAX_DIMENSIONS];
      long_int GlobalLeftEdgeOther[MAX_DIMENSIONS];
      int dim;
      this->GetLeftEdge(GlobalLeftEdgeThis);
      Other->GetLeftEdge(GlobalLeftEdgeOther);
      Other->GetCellDimensions(OtherDims);
      for (dim = 0; dim < this->Rank; dim++) {
        OverlapStart = MAX(GlobalLeftEdgeThis[dim], GlobalLeftEdgeOther[dim]);
        OverlapEnd = MIN(GlobalLeftEdgeThis[dim] + this->CellDimensions[dim],
                         GlobalLeftEdgeOther[dim] + OtherDims[dim]);
        LeftEdgeThis[dim] = OverlapStart - GlobalLeftEdgeThis[dim];
        LeftEdgeOther[dim] = OverlapStart - GlobalLeftEdgeOther[dim];
        CopyDims[dim] = MAX(OverlapEnd - OverlapStart, 0);
      }
}

void FieldDescriptor::PrintFieldInformation() {
  fprintf(stdout,"\nField: %s\n",this->Name);
  fprintf(stdout,"  UnitsName: %s\n",this->UnitsName);
  fprintf(stdout,"  CellDimensions: ");
  for (int dim=0; dim<MAX_DIMENSIONS; dim++)
    fprintf(stdout," %"ISYM"", this->CellDimensions[dim]);
  fprintf(stdout,"\n  LeftEdge: ");
  for (int dim=0; dim<MAX_DIMENSIONS; dim++)
    fprintf(stdout," %li", this->LeftEdge[dim]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  Level: %"ISYM"\n",this->Level);
  fprintf(stdout,"  DeallocateFieldPointer: %"ISYM"\n",this->DeallocateFieldPointer);
  fprintf(stdout,"  DeallocateFieldvalues: %"ISYM"\n",this->DeallocateFieldValues);
}


// Mathematical Operations

float FieldDescriptor::Min() {
  return this->Min(NULL, NULL);
}

float FieldDescriptor::Min(int *LeftEdge, int *RightEdge) {
  return this->UnaryAccumulator<MinVal>(
      LeftEdge, RightEdge, 1e300);
}

float FieldDescriptor::Max() {
  return this->Max(NULL, NULL);
}

float FieldDescriptor::Max(int *LeftEdge, int *RightEdge) {
  return this->UnaryAccumulator<MaxVal>(
      LeftEdge, RightEdge, -1e300);
}

float FieldDescriptor::Sum() {
  return this->Sum(NULL, NULL);
}

float FieldDescriptor::Sum(int *LeftEdge, int *RightEdge) {
  float v = this->UnaryAccumulator<AddVal>(
      LeftEdge, RightEdge, 0.0);
}

// Operations from other FieldDescriptors

void FieldDescriptor::CopyFrom(FieldDescriptor *Other) {
  int LeftEdgeThis[MAX_DIMENSIONS];
  int LeftEdgeOther[MAX_DIMENSIONS];
  int CopyDims[MAX_DIMENSIONS];
  // This only considers *same-level* copies
  this->GetOverlapRegion(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
  this->InPlaceBinaryOperation<CopyVal>(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
}

void FieldDescriptor::CopyFrom(float val) {
  static int Zero[MAX_DIMENSIONS] = {0, 0, 0};
  this->InPlaceBinaryOperation<CopyVal>(val, Zero, this->CellDimensions);
}

void FieldDescriptor::Add(FieldDescriptor *Other) {
  int LeftEdgeThis[MAX_DIMENSIONS];
  int LeftEdgeOther[MAX_DIMENSIONS];
  int CopyDims[MAX_DIMENSIONS];
  // This only considers *same-level* copies
  this->GetOverlapRegion(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
  this->InPlaceBinaryOperation<AddVal>(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
}

void FieldDescriptor::Add(float val) {
  static int Zero[MAX_DIMENSIONS] = {0, 0, 0};
  this->InPlaceBinaryOperation<AddVal>(val, Zero, this->CellDimensions);
}

void FieldDescriptor::Subtract(FieldDescriptor *Other) {
  int LeftEdgeThis[MAX_DIMENSIONS];
  int LeftEdgeOther[MAX_DIMENSIONS];
  int CopyDims[MAX_DIMENSIONS];
  // This only considers *same-level* copies
  this->GetOverlapRegion(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
  this->InPlaceBinaryOperation<SubVal>(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
}

void FieldDescriptor::Subtract(float val) {
  static int Zero[MAX_DIMENSIONS] = {0, 0, 0};
  this->InPlaceBinaryOperation<SubVal>(val, Zero, this->CellDimensions);
}

void FieldDescriptor::Multiply(FieldDescriptor *Other) {
  int LeftEdgeThis[MAX_DIMENSIONS];
  int LeftEdgeOther[MAX_DIMENSIONS];
  int CopyDims[MAX_DIMENSIONS];
  // This only considers *same-level* copies
  this->GetOverlapRegion(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
  this->InPlaceBinaryOperation<MultVal>(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
}

void FieldDescriptor::Multiply(float val) {
  static int Zero[MAX_DIMENSIONS] = {0, 0, 0};
  this->InPlaceBinaryOperation<MultVal>(val, Zero, this->CellDimensions);
}

void FieldDescriptor::Divide(FieldDescriptor *Other) {
  int LeftEdgeThis[MAX_DIMENSIONS];
  int LeftEdgeOther[MAX_DIMENSIONS];
  int CopyDims[MAX_DIMENSIONS];
  // This only considers *same-level* copies
  this->GetOverlapRegion(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
  this->InPlaceBinaryOperation<DivVal>(Other, LeftEdgeThis, LeftEdgeOther, CopyDims);
}

void FieldDescriptor::Divide(float val) {
  static int Zero[MAX_DIMENSIONS] = {0, 0, 0};
  this->InPlaceBinaryOperation<DivVal>(val, Zero, this->CellDimensions);
}

// Operator Overloading

FieldDescriptor &FieldDescriptor::operator+=(FieldDescriptor *Other) {
    this->Add(Other);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator+=(float val) {
    this->Add(val);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator-=(FieldDescriptor *Other) {
    this->Subtract(Other);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator-=(float val) {
    this->Subtract(val);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator*=(FieldDescriptor *Other) {
    this->Multiply(Other);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator*=(float val) {
    this->Multiply(val);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator/=(FieldDescriptor *Other) {
    this->Divide(Other);
    return (*this);
}

FieldDescriptor &FieldDescriptor::operator/=(float val) {
    this->Divide(val);
    return (*this);
}

// Centering Stuff

int FieldDescriptor::CanCombine(FieldDescriptor *Other) {
    if (Other->ValueCentering != this->ValueCentering) return 0;
    if (Other->Rank != this->Rank) return 0;
    return 1;
}

// Operations

template <MathFunction function>
float FieldDescriptor::UnaryAccumulator(
    int *LeftEdge, int *RightEdge, float InitialValue) {
  // NOTE: This takes a RightEdge in *cell* values.  This will *do the right
  // thing* for face and vertex centered fields.
  int i, j, k;
  float val = InitialValue;
  int ind;
  static int Zero[MAX_DIMENSIONS] = {0, 0, 0};
  int Extension[MAX_DIMENSIONS];
  if (LeftEdge == NULL) {
    LeftEdge = Zero;
  }
  if (RightEdge == NULL) {
    RightEdge = this->CellDimensions;
  }
  this->GetFieldExtension(Extension);
  float *vals = this->GetValues();

  for (i = LeftEdge[0]; i < RightEdge[0] + Extension[0]; i++) {
    for (j = LeftEdge[1]; j < RightEdge[1] + Extension[1]; j++) {
      for (k = LeftEdge[2]; k < RightEdge[2] + Extension[2]; k++) {
        ind = this->Index(i, j, k);
        val = function(val, vals[ind]);
      }
    }
  }
  return val;
};

template <MathFunction function>
void FieldDescriptor::InPlaceBinaryOperation(
    FieldDescriptor *Other,
    int LeftEdgeThis[MAX_DIMENSIONS],
    int LeftEdgeOther[MAX_DIMENSIONS],
    int CopyDims[MAX_DIMENSIONS]) {

  if (this->CanCombine(Other) == 0) {
    throw FieldsIncompatible("FieldDescriptors are not compatible.");
  }

  int Extension[MAX_DIMENSIONS];
  this->GetFieldExtension(Extension);

  int dim, i, j, k, i1, i2, j1, j2, k1, k2, ind1, ind2;

  float *v1 = this->GetValues();
  float *v2 = Other->GetValues();

  for (i = 0, i1 = LeftEdgeThis[0], i2 = LeftEdgeOther[0];
       i < CopyDims[0] + Extension[0];
       i++, i1++, i2++) {
    for (j = 0, j1 = LeftEdgeThis[1], j2 = LeftEdgeOther[1];
        j < CopyDims[1] + Extension[1];
        j++, j1++, j2++) {
      for (k = 0, k1 = LeftEdgeThis[2], k2 = LeftEdgeOther[2];
          k < CopyDims[2] + Extension[2];
          k++, k1++, k2++) {
        ind1 = this->Index(i1, j1, k1);
        ind2 = Other->Index(i2, j2, k2);
        v1[ind1] = function(v1[ind1], v2[ind2]);
      }
    }
  }
}

template <MathFunction function>
void FieldDescriptor::InPlaceBinaryOperation(
    float OtherValue,
    int LeftEdge[MAX_DIMENSIONS],
    int CopyDims[MAX_DIMENSIONS]) {

  int dim, i, j, k, i1, j1, k1, ind1;

  float *v1 = this->GetValues();

  int Extension[MAX_DIMENSIONS];
  this->GetFieldExtension(Extension);

  for (i = 0, i1 = LeftEdge[0]; i < CopyDims[0] + Extension[0]; i++, i1++) {
    for (j = 0, j1 = LeftEdge[1]; j < CopyDims[1] + Extension[1]; j++, j1++) {
      for (k = 0, k1 = LeftEdge[2]; k < CopyDims[2] + Extension[2]; k++, k1++) {
        ind1 = this->Index(i1, j1, k1);
        v1[ind1] = function(v1[ind1], OtherValue);
      }
    }
  }
}
