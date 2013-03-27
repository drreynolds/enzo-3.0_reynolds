//
// Field Registry
//
// Authors: Matthew Turk
//          Greg Bryan

#ifdef min
#define OLDmin min
#undef min
#endif

#ifdef max
#define OLDmax max
#undef max
#endif

#include <map>

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

#ifdef OLDmin
#define min OLDmin
#endif

#ifdef OLDmax
#define max OLDmax
#endif
