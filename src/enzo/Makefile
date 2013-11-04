#=======================================================================
#
# FILE:        Makefile.config
#
# SUMMARY:     Configurable Makefile for Enzo
#
# AUTHOR:      James Bordner (jobordner@ucsd.edu)
#
# DATE:        2007-02-21
#
# DESCRIPTION
#              See 'gmake help' for definitive description of targets
#
#              Makefile.config includes the following files:
# 
#              Make.config.settings   default configuration settings
#              Make.config.override   local user settings
#              Make.config.assemble   maps 'config' settings to 'flag' settings
#              Make.config.targets    configure targets
#              Make.mach.*            all machine-dependent settings
#              Make.config.objects    list of object files
#              DEPEND                 Make-generated dependencies
#
#              Make.mach.* should be the only file that one should
#              edit when porting Enzo to a new machine.
#
#              Make.config.override should be the only file that
#              one should edit when defining local user settings.
#              Preferably, this is done implicitly through
#              the available make targets (e.g. "gmake precision-32").
#              Use "gmake help-config" for a list of all configuration
#              settings.  These make targets do error-checking; hand-editing 
#              Make.config.override is more error-prone.
#
#=======================================================================

# Use bash since sh on datastar does not recognize ">&" used in dep: target

SHELL    = /bin/bash

TOP_DIR  = ../..
EXE      = enzo
OUTPUT   = out.compile

ENZO_DIR = .
CONFIG_DIR = build/config
BUILD_DIR = build/objects
GEN_DIR = build/generated_code
UTILS_DIR = build/utilities
MODULES  = 
VERBOSE  = 0
ETAGS_COMMAND = "0"

LFLAGS=--header-file=libconfig/scanner.h --prefix=libconfig_yy

SVN      = hg

vpath %.o $(BUILD_DIR)/
vpath %.C $(GEN_DIR)/
vpath %.h $(GEN_DIR)/

#-----------------------------------------------------------------------
# Make.config.settings is used for setting default values to all compile-time 
# configuration settings.
#-----------------------------------------------------------------------

include $(ENZO_DIR)/$(CONFIG_DIR)/Make.config.settings

#-----------------------------------------------------------------------
# Make.config.machine is used for setting which Make.mach.* file to use
#-----------------------------------------------------------------------

MAKE_CONFIG_MACHINE  = $(ENZO_DIR)/$(CONFIG_DIR)/Make.config.machine
include $(ENZO_DIR)/$(CONFIG_DIR)/Make.config.machine

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing Enzo with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(ENZO_DIR)/$(CONFIG_DIR)/Make.config.override

include $(MAKE_CONFIG_OVERRIDE)

# THIS WAY OF DOING THE ABOVE DOES NOT WORK:
#
# include $(ENZO_DIR)/$(CONFIG_DIR)/Make.config.override

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include $(CONFIG_DIR)/Make.config.assemble

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#-----------------------------------------------------------------------

MACH_SHARED_EXT=so

-include $(ENZO_DIR)/build/machines/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.enzo/Make.mach.$(CONFIG_MACHINE)

#=======================================================================
# OBJECT FILES
#=======================================================================

include io/python_bridge/Make.module.objects
include communication/Make.module.objects
include control/Make.module.objects
include grid/Make.module.objects
include hierarchy/Make.module.objects
include hydro_rk/Make.module.objects
include io/Make.module.objects
include particles/Make.module.objects
include problem_types/Make.module.objects
include solvers/Make.module.objects
include utilities/Make.module.objects

OBJS_CONFIG_LIB=$(ENZO_OBJS)
OBJS_LINK=$(addprefix $(BUILD_DIR)/, $(OBJS_LIB))

#-----------------------------------------------------------------------
# MAKE ENZO BY DEFAULT
#-----------------------------------------------------------------------

all: MKDIR MACHNOTES MAKETAGS $(EXE).exe
.DEFAULT_GOAL: all

# I believe this has to be done here, as it is target-specific.
lib: CFLAGS   += $(MACH_SHARED_FLAGS) 
lib: CXXFLAGS += $(MACH_SHARED_FLAGS) -DSHARED_LIBRARY
lib: FFLAGS   += $(MACH_SHARED_FLAGS)
lib: F90FLAGS += $(MACH_SHARED_FLAGS)
lib: LDFLAGS  += $(MACH_SHARED_FLAGS)
lib: lib$(EXE)_p$(ASSEMBLE_PARTICLE_NUMBER)_b$(ASSEMBLE_PRECISION_NUMBER).$(MACH_SHARED_EXT)

#-----------------------------------------------------------------------
# MAKE AN EXECUTABLE
#-----------------------------------------------------------------------

.PHONY: MKDIR
MKDIR: $(BUILD_DIR) $(GEN_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(GEN_DIR):
	mkdir -p $(GEN_DIR)

MAKETAGS:
	@(if [ "$(ETAGS_COMMAND)" != "0" ]; then \
	   	echo -e "Making Tags with $(ETAGS_COMMAND)" ; \
		rm -f TAGS ; \
		$(ETAGS_COMMAND) ; \
		echo -e "Done updating TAGS" ; \
	fi)
MACHNOTES: 
	@echo -e $(MACHINE_NOTES)

%.exe: $(MODULES) autogen dep %.o $(OBJS_LIB) MACHNOTES
	@rm -f $@
	@echo "Linking enzo executable. Type  cat $(OUTPUT)  in case it fails."
	@(if [ $(VERBOSE) -eq 0 ]; then \
		$(LD) $(LDFLAGS) -o $*.exe $(BUILD_DIR)/$*.o $(OBJS_LINK) $(LIBS) >& $(OUTPUT) ; \
	else \
		$(LD) $(LDFLAGS) -o $*.exe $(BUILD_DIR)/$*.o $(OBJS_LINK) $(LIBS) >> $(OUTPUT) \
		2>&1 ; \
	fi)	
	@(if [ -e $@ ] ; then  \
		echo "Success!" ; \
		echo "Compiled enzo from " ; \
		echo "Mercurial Branch   `$(SVN) identify -b`" ; \
		echo "Mercurial Revision `$(SVN) identify -i`" ; \
		if [ ! -e $(TOP_DIR)/bin ]; then mkdir $(TOP_DIR)/bin; fi; \
		if [ "$(CONFIG_OPT)" == "debug" ]; then \
		echo "Compiled with opt-debug (OPT_FLAGS = $(ASSEMBLE_OPT_FLAGS))."; \
		echo "May run much faster while still accurate with make opt-high (OPT_FLAGS = $(MACH_OPT_HIGH))"; \
		fi; \
                cp $(EXE).exe $(TOP_DIR)/bin/$(EXE); \
	else \
		echo "$(LD) $(LDFLAGS) -o $*.exe $(BUILD_DIR)/$*.o $(OBJS_LINK) $(LIBS)" > temp1; \
		cat temp1 $(OUTPUT) > temp2 ;\
		mv -f temp2 $(OUTPUT); \
		rm -f temp1 temp2; \
		echo "Failed! See $(OUTPUT) for error messages"; \
	fi)

lib%.$(MACH_SHARED_EXT): $(MODULES) autogen dep enzo.o $(OBJS_LIB)
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) $(SHARED_OPT) -o $@ enzo.o $(OBJS_LINK) $(LIBS) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
	   echo "Success!"; \
	else \
	   echo "$(LD) $(LDFLAGS) $(SHARED_OPT) -o $*.so $(BUILD_DIR)/$*.o $(OBJS_LINK) $(LIBS)" >> temp1; \
	   cat temp1 $(OUTPUT) > temp2; \
	   rm -f temp1; \
	   mv -f temp2 $(OUTPUT); \
	   echo "Failed! See $(OUTPUT) for error messages"; \
	fi)

io/InitializePythonInterface.C: io/InitializePythonInterface_finderfunctions.inc
io/InitializePythonInterface_finderfunctions.inc: headers/global_data.h
	@python $(UTILS_DIR)/create_dictionary_mapping.py 

#-----------------------------------------------------------------------
# WRITE ALL COMPILER OUTPUT TO FILE
#-----------------------------------------------------------------------

.PHONY: verbose
verbose: VERBOSE = 1
#verbose:
#	@rm -fv $(OUTPUT)
#	@echo "Writing all compiler output to $(OUTPUT)"
verbose: $(EXE).exe

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .F .F90 .o .cu .pyx

# Inhibit removing any *.o files after compiling

.PRECIOUS: %.o %.C

#.src.f:
#	@$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $<
.F.o:
	@echo "Compiling $<"
	@rm -f $(BUILD_DIR)/$@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(FC) -c -o $(BUILD_DIR)/$@ $(FFLAGS) -DFILENAME=$*.F $(DEFINES) $(INCLUDES) $< >& $(OUTPUT) ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
             echo; \
             echo "$(FC) -c -o $(BUILD_DIR)/$@ $(FFLAGS) -DFILENAME=$*.F $(DEFINES) $(INCLUDES) $<"; \
             echo; \
             $(FC) -c -o $(BUILD_DIR)/$@ $(FFLAGS) -DFILENAME=$*.F $(DEFINES) $(INCLUDES) $<; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(FC) -c -o $(BUILD_DIR)/$@ $(FFLAGS) -DFILENAME=$*.F $(DEFINES) $(INCLUDES) $< >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

#.src90.f90:
#	@echo "Compiling $<"
#	@$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $<
.F90.o:
	@rm -f $(BUILD_DIR)/$@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(F90) -c -o $(BUILD_DIR)/$@ $(F90FLAGS) -DFILENAME=$*.F90 $(DEFINES) $(INCLUDES) $< >& $(OUTPUT) ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
             echo; \
             echo "$(F90) -c -o $(BUILD_DIR)/$@ -DFILENAME=$*.F90 $(DEFINES) $(F90FLAGS) $(INCLUDES) $<"; \
             echo; \
             $(F90) -c -o $(BUILD_DIR)/$@ -DFILENAME=$*.F90 $(DEFINES) $(F90FLAGS) $(INCLUDES) $<; \
             echo; \
             exit 1; \
	  fi ; \
	else \
	  $(F90) -c -o $(BUILD_DIR)/$@ $(F90FLAGS) -DFILENAME=$*.F90 $(DEFINES) $(INCLUDES) $< >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.c.o:
	@rm -f $(BUILD_DIR)/$@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CC) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CFLAGS) $(INCLUDES) $< \
	    >& $(OUTPUT) ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
             echo; \
             echo "$(CC) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CFLAGS) $(INCLUDES) $<"; \
             echo; \
             $(CC) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CFLAGS) $(INCLUDES) $<;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CC) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CFLAGS) $(INCLUDES) $< \
	    >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.pyx.o:
	@rm -f $(BUILD_DIR)/$@
	@echo "Compiling $<"
	@cython --cplus $< -o $(GEN_DIR)/$*.C
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $(GEN_DIR)/$*.C \
	    >& $(OUTPUT) ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
             echo; \
             echo "$(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C"; \
             echo; \
             $(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $<;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $(GEN_DIR)/$*.C \
	    >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.C.o:
	@rm -f $(BUILD_DIR)/$@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $< \
	    >& $(OUTPUT) ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
             echo; \
             echo "$(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $<"; \
             echo; \
             $(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $<;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CXX) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $< \
	    >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.cu.o: 
	@rm -f $(BUILD_DIR)/$@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CUDACOMPILER) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(INCLUDES) $(CUDACOMPFLAGS) \
	    $(INCLUDES) $< >& $(OUTPUT) ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo; \
             echo "$(CUDACOMPILER) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CUDACOMPFLAGS) $(INCLUDES) $<"; \
             echo; \
             $(CUDACOMPILER)  -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CUDACOMPFLAGS) $(INCLUDES) $<;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CUDACOMPILER) -c -o $(BUILD_DIR)/$@ $(DEFINES) $(CUDACOMPFLAGS) \
	    $(INCLUDES) $< >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $(BUILD_DIR)/$@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

#-----------------------------------------------------------------------
# Generate all make-generated source files
#-----------------------------------------------------------------------

.PHONY: autogen
autogen: io/auto_show_config.C io/auto_show_flags.C io/auto_show_version.C io/auto_show_compile_options.C

# Force update of auto_show_config.C

.PHONY: io/auto_show_config.C
io/auto_show_config.C:
	-@$(MAKE) -s show-config  >& temp.show-config
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_config(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-config > io/auto_show_config.C

# Force update of auto_show_flags.C

.PHONY: io/auto_show_flags.C
io/auto_show_flags.C:
	-@$(MAKE) -s show-flags  >& temp.show-flags
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_flags(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-flags > io/auto_show_flags.C

# Force update of auto_show_version.C

.PHONY: io/auto_show_version.C
io/auto_show_version.C:
	-@$(MAKE) -s show-version  >& temp.show-version
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_version(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-version > io/auto_show_version.C

# Force update of auto_show_compile_options.C

.PHONY: io/auto_show_compile_options.C
io/auto_show_compile_options.C:
	-@python $(UTILS_DIR)/create_config_info.py

#-----------------------------------------------------------------------
# Generate dependency file
#-----------------------------------------------------------------------

.PHONY: dep
dep:
	@echo "Updating DEPEND"
	-@(makedepend $(DEFINES) $(INCLUDES) -fDEPEND -o.o -m -- -- */*.C  ) >& out.make.DEPEND	
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.C) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.c) >> out.make.DEPEND 2>&1	
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.F) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.F90) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- */*.h) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.h) >> out.make.DEPEND 2>&1

include DEPEND

#-----------------------------------------------------------------------
# Radiative transfer module
#-----------------------------------------------------------------------

#include $(ENZO_DIR)/photons/$(CONFIG_DIR)/Make.config.objects
#
#.PHONY: photon
#photon: OBJS_LIB += photons/*.o
#photon:
#	@echo "Making radiative transfer module"
#	+(cd photons/ ; make photon)

#-----------------------------------------------------------------------
# HELP TARGET
#-----------------------------------------------------------------------

help:
	@echo
	@echo "========================================================================"
	@echo "   Enzo Makefile Help"
	@echo "========================================================================"
	@echo
	@echo "   gmake                Compile and generate the executable 'enzo.exe'"
	@echo "   gmake install        Copy the executable to bin/enzo"
	@echo "   gmake help           Display this help information"
	@echo "   gmake clean          Remove object files, executable, etc."
	@echo "   gmake dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   gmake show-version   Display revision control system branch and revision"
	@echo "   gmake show-diff      Display local file modifications"
	@echo
	@echo "   gmake help-config    Display detailed help on configuration make targets"
	@echo "   gmake show-config    Display the configuration settings"
	@echo "   gmake show-flags     Display specific compilation flags"
	@echo "   gmake default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------

clean:
	-@rm -f *.so *.o uuid/*.o *.mod *.f *.f90 DEPEND.bak *~ $(OUTPUT) *.exe \
          $(BUILD_DIR)/*.o \
          $(BUILD_DIR)/*.C \
          $(BUILD_DIR)/*.h \
          $(GEN_DIR)/*.C $(GEN_DIR)/*.h \
          io/auto_show*.C hydro_rk/*.o *.oo hydro_rk/*.oo \
          uuid/*.oo DEPEND TAGS \
          libconfig/*.o 
	-@touch DEPEND

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(ENZO_DIR)/$(CONFIG_DIR)/Make.config.targets
# DO NOT DELETE
