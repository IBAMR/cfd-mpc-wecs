######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
IBAMR_SRC_DIR   = /path/to/IBAMR_SRC_DIR
IBAMR_BUILD_DIR = /path/to/IBAMR_BUILD_DIR

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

## Needed for Xcode to capture compiler errors and warnings.
ifdef XCODE_VERSION_ACTUAL
CXXFLAGS += -fno-color-diagnostics
endif

######################################################################
## Build the application.
##
## NOTE: The following assumes that all .cpp files in the present
##       directory are used to build the executable.
SRC = $(wildcard *.cpp)
CPPFLAGS += -MD -MP
PDIM = 3

OBJS3D = $(SRC:%.cpp=%.o) $(IBAMR_LIB_3D) $(IBTK_LIB_3D)
main3d: $(OBJS3D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o $@
clean:
	$(RM) main3d
	$(RM) *.o *.lo *.objs *.ii *.int.c *.d
	$(RM) -r .libs

-include $(SRC:%.cpp=%.d)

