"""
The necessary imports for accessing Enzo data

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://enzo-project.org/
License: This file is covered under the Enzo license
"""

ctypedef double Eflt
ctypedef double FLOAT

# Now the business with the ints
ctypedef long long long_int
ctypedef long_int Eint
ctypedef int Eint32
ctypedef long_int Eint64
ctypedef long_int PINT

cdef extern from "preincludes.h":
    pass

cdef extern from "math.h":
    pass

cdef extern from "string.h":
    pass

cdef extern from "map":
    pass

cdef extern from "stdio.h":
    pass

cdef extern from "math.h":
    pass

cdef extern from "iostream":
    pass

cdef extern from "cmath":
    pass

cdef extern from "complex":
    pass

cdef extern from "communicators.h":
    pass

cdef extern from "performance.h":
    pass

cdef extern from "macros_and_parameters.h":
    pass

cdef extern from "typedefs.h":
    pass

cdef extern from "global_data.h":
    pass

cdef extern from "Fluxes.h":
    pass

cdef extern from "GridList.h":
    pass

cdef extern from "ExternalBoundary.h":
    pass

cdef extern from "Grid.h":
    cdef cppclass grid:
        pass

cdef extern from "Hierarchy.h":
    pass

cdef extern from "communication.h":
    pass

cdef extern from "CommunicationUtilities.h":
    pass


