"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://enzo-project.org/
License: This file is covered under the Enzo license
"""

import numpy as np
cimport numpy as np
cimport cython
cimport libc.stdlib
from libcpp.string cimport string
from cython.operator cimport dereference as deref

from enzo_includes cimport *

cdef extern from "EventHooks.h":
    void RunEventHooks(string,
                       c_HierarchyEntry **Grids,
                       c_TopGridData &MetaData,
                       EventDataContainer *LocalData)

cdef extern from "EventDataContainers.h":
    cdef cppclass EventDataContainer:
        pass

    cdef cppclass TestEventDataContainer:
        int SomethingToPrint

cdef class MessageCoordinator:
    cdef c_HierarchyEntry **Grids
    cdef c_TopGridData MetaData

    def issue_test_event(self, value_to_print):
        cdef TestEventDataContainer *pedc = new TestEventDataContainer()
        pedc.SomethingToPrint = value_to_print
        cdef EventDataContainer *data = <EventDataContainer *> pedc
        RunEventHooks("Python Event".encode("UTF-8"), self.Grids, self.MetaData, data)
        del pedc

cdef MessageCoordinator mc = MessageCoordinator()
test_event = mc.issue_test_event

cdef public void create_coordinator(c_HierarchyEntry **Grids,
                                    c_TopGridData &MetaData):
    mc.Grids = Grids
    mc.MetaData = MetaData

cdef extern from "fix_enzo_defs.h":
    pass
