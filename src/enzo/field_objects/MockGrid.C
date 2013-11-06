//
// Field Descriptor Object
//
// Authors: Matthew Turk
//          Greg Bryan

#include "FieldObjects.h"
#include "malloc.h"
#include "assert.h"

Grid::Grid() {
  this->NumberOfFields = 0;
  this->Dimensions[0] = 16;
  this->Dimensions[1] = 16;
  this->Dimensions[2] = 16;
  this->LeftEdge[0] = 0;
  this->LeftEdge[1] = 0;
  this->LeftEdge[2] = 0;

  int field;
  for (field = 0; field < MAX_FIELDS; field++) {
    this->BaryonFields[field] = NULL;
  }
}

Grid::Grid(int Dimensions[MAX_DIMENSIONS],
           long_int LeftEdge[MAX_DIMENSIONS]) {
  int field;
  this->NumberOfFields = 0;
  this->Dimensions[0] = Dimensions[0];
  this->Dimensions[1] = Dimensions[1];
  this->Dimensions[2] = Dimensions[2];
  this->LeftEdge[0] = LeftEdge[0];
  this->LeftEdge[1] = LeftEdge[1];
  this->LeftEdge[2] = LeftEdge[2];
  for (field = 0; field < MAX_FIELDS; field++) {
    this->BaryonFields[field] = NULL;
  }
}

Grid::~Grid() {
  FieldRegistry::iterator iter;
  for (iter = this->Fields.begin();
       iter != this->Fields.end(); ++iter) {
    if (iter->second != NULL) {
      delete (iter->second);
    }
  }
}

FieldDescriptor* Grid::operator[](const char *Name) {
  return this->GetField(Name);
};

FieldDescriptor* Grid::GetField(const char *Name) {
  FieldDescriptor *fd = this->Fields[Name];
  assert(fd != NULL);
  return fd;
};

void Grid::AttachField(char * Name, FieldDescriptor *NewField) {
  // This does not get added to BaryonFields
  assert(NewField != NULL);
  this->Fields[Name] = NewField;
}

FieldDescriptor* Grid::AddField(char *Name, char *UnitsName, int NumberOfGhostZones,
                          InterpolationType InterpolationMethod,
                          CenteringType ValueCentering) {
  // We own this field.
  int Dimensions[MAX_DIMENSIONS];
  long_int LeftEdge[MAX_DIMENSIONS];
  int dim;
  for (dim = 0; dim < MAX_DIMENSIONS; dim++) {
    Dimensions[dim] = this->Dimensions[dim] + 2*NumberOfGhostZones;
    LeftEdge[dim] = this->LeftEdge[dim] - NumberOfGhostZones;
  }
  FieldDescriptor *NewField = new FieldDescriptor(
        ValueCentering, 3, Dimensions, LeftEdge,
        InterpolationMethod, Name, UnitsName,
        this->BaryonFields + this->NumberOfFields++);
  this->Fields[Name] = NewField;
  return NewField;
}


int Grid::FindField(int FieldType) {
    for (int i = 0; i < NumberOfFields; i++) {
      if (FieldType == this->FieldTypes[i]) return i;
    }
}
