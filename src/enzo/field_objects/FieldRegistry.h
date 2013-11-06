//
// Field Registry
//
// Authors: Matthew Turk
//          Greg Bryan
#ifndef __FIELD_REGISTRY_H__
#define __FIELD_REGISTRY_H__

typedef std::map<std::string, FieldDescriptor *> FieldRegistry;
typedef std::map<int, std::string> FieldNumbers;
void FillFieldRegistry(int Rank, FieldRegistry *Fields, FieldNumbers *FieldIDs);

#ifdef DEFINE_STORAGE
# define EXTERN_FIELDS
#else /* DEFINE_STORAGE */
# define EXTERN_FIELDS extern
#endif

EXTERN_FIELDS FieldRegistry BaseFieldTypes;
EXTERN_FIELDS FieldNumbers BaseFieldIDs;

#undef EXTERN_FIELDS
#endif
