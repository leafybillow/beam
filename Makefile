ROOTCONFIG = root-config

CXX           = $(shell $(ROOTCONFIG) --cxx)
# CXXFLAGS      = -O0 -pipe -Wall -W -Woverloaded-virtual -Wextra -Wno-missing-field-initializers -fPIC -ldl
CXXFLAGS	= -O0 -Wall -Woverloaded-virtual -fPIC -Wextra
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)
ROOTINC :=$(shell $(ROOTCONFIG) --incdir)
CXXFLAGS+= $(ROOTCFLAGS)
CXXFLAGS+= -g
LIBS	:= $(ROOTLIBS)
LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LD	:= $(shell $(ROOTCONFIG) --ld)
OBJS	:= BeamGEMPlane.o BeamGEMTracker.o
HDR	:= $(OBJS:.o=.h)

all:  $(OBJS) beam_Dict beam

beam_Dict: $(HDR) beam_LinkDef.h 
	rootcint -f $@.cc -c -I$(ROOTINC) $^;
	$(CXX) $(CXXFLAGS) -c -o $@.o $@.cc ;
$(OBJS):
	$(CXX) $(CXXFLAGS) -c -o $@ $(@:.o=.cc) ;

beam: $(OBJS) beam_Dict.o
	$(LD) $(CXXFLAGS) $(LDFLAGS) -shared $^ -o lib$@.so;
clean:
	rm *.o;
	rm *Dict*;	
	rm *.so;	