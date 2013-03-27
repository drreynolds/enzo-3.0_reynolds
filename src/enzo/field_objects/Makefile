GTEST_DIR=gtest
CC=g++
DEFINES=-DGTEST_HAS_PTHREAD=0 -DMOCK_GRID
CCFLAGS=-ggdb -O0 -Ifield_objects/ -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(DEFINES)
LDFLAGS=-L$(GTEST_DIR)/lib 
TEST_SRC=TestFieldDescriptors.C \
		 TestGrid.C \
		 TestFieldCenterings.C \
		 TestFaceCenteredFields.C \
		 TestEdgeCenteredFields.C \
		 TestCornerCenteredFields.C \
		 TestFieldBases.C
TEST_OBJ=$(TEST_SRC:.C=.o)
SRC=FieldDescriptor.C \
	MockGrid.C \
	FieldRegistry.C
OBJ=$(SRC:.C=.o)
EXE=TestRunner.exe

%.exe: %.o $(OBJ) $(TEST_OBJ)
	    g++ $(CCFLAGS) $(LDFLAGS) $(TEST_OBJ) $(OBJ) ${GTEST_DIR}/src/gtest-all.cc $< -o $@
%.o: %.C
	    g++ $(CCFLAGS) -c $< -o $@

all: $(OBJ) $(EXE)

.PHONY: clean

runtests: $(TEST_OBJ) $(EXE)
	./$(EXE)

clean:
	rm -f $(TEST_OBJ)
	rm -f $(OBJ)
	rm -f $(EXE)
