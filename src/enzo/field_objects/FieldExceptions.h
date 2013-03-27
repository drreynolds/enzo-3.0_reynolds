
//
// Field Object Exceptions
//
// Authors: Matthew Turk
//          Greg Bryan

#include <stdexcept>
#include <string>

class FieldsIncompatible : public std::runtime_error {
public:
    FieldsIncompatible(const std::string& message) 
        : std::runtime_error(message) { };
};
