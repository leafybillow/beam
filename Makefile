ROOTCONFIG = root-config

CXX     =$(shell $(ROOTCONFIG) --cxx)
GCC	=$(shell $(ROOTCONFIG) --cc)
CXXFLAGS	= -O0 -Wall -Woverloaded-virtual -fPIC -Wextra
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs) -lMinuit
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)
ROOTINC :=$(shell $(ROOTCONFIG) --incdir)

MAKEDEPEND =$(CXX)
INCLUDES = -I$(ROOTINC)

CXXFLAGS+= $(ROOTCFLAGS)
CXXFLAGS+= -g

LIBS	:= $(ROOTLIBS)
LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LD	:= $(shell $(ROOTCONFIG) --ld)
OBJS	:= \
	BeamGEMPlane.o BeamGEMTracker.o BeamGEMStrip.o\
	BeamGEMProjection.o BeamGEMData.o\
	BeamConfig.o BeamAnalysis.o

HDR	:= $(OBJS:.o=.h) BeamParameters.h BeamTypes.h
DEPS 	:= $(OBJS:.o=.d)

all:  $(OBJS) $(BINARIES) beam_Dict libbeam beam

beam_Dict: $(HDR) beam_LinkDef.h
	rootcint -f $@.cc -c -I$(ROOTINC) $^;
	$(CXX) $(CXXFLAGS) -c -o $@.o $@.cc ;
$(OBJS):
	$(CXX) $(CXXFLAGS) -c -o $@ $(@:.o=.cc) ;

libbeam: $(OBJS) beam_Dict.o
	$(LD) $(CXXFLAGS) $(LDFLAGS) -shared $^ -o $@.so;

beam:	$(OBJS) beam_Dict.o
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $^ $@.cc -o $@;
clean:
	rm -f *.o;
	rm -f *Dict*;	
	rm -f *.so;	

%.d:	%.cc
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -MM $(INCLUDES) -c $< \
	      | sed '\''s%^.*\.o%$*\.o%g'\'' \
	      | sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'

-include $(DEPS)
