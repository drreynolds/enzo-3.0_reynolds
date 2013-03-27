//
// Field Registry
//
// Authors: Matthew Turk
//          Greg Bryan

#include <map>

typedef std::map<std::string, FieldDescriptor *> FieldRegistry;
void FillFieldRegistry(int Rank, FieldRegistry *Fields);

