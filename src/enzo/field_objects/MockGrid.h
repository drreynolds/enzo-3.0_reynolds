//
// Mock Objects needed for testing
//
// Authors: Matthew Turk
//          Greg Bryan

#ifndef __MOCK_OBJECTS_H__
#define __MOCK_OBJECTS_H__

#include <string>
#include <map>

class Grid
{
    public:
      Grid();
      Grid(int Dimensions[MAX_DIMENSIONS],
           long_int LeftEdge[MAX_DIMENSIONS]);
      ~Grid();

      FieldDescriptor *operator[](const char *Name);
      FieldDescriptor *GetField(const char *Name);

      void AttachField(char* Name, FieldDescriptor *NewField);
      FieldDescriptor *AddField(char* Name, char *UnitsName,
               int NumberOfGhostZones,
               InterpolationType InterpolationMethod,
               CenteringType ValueCentering);

      float *BaryonFields[MAX_FIELDS];
      int NumberOfFields;
      int Dimensions[MAX_DIMENSIONS];
      long_int LeftEdge[MAX_DIMENSIONS];

      FieldRegistry Fields;
      int FindField(int FieldType);
      int FieldTypes[MAX_FIELDS];

    friend class FieldCollection;

    protected:

    private:

};

#endif
