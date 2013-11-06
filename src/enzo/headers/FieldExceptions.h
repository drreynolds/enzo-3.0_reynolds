
//
// Field Object Exceptions
//
// Authors: Matthew Turk
//          Greg Bryan

class FieldsIncompatible : public std::runtime_error {
public:
    FieldsIncompatible(const std::string& message) 
        : std::runtime_error(message) { };
};
