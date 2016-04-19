Bindir = .

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

ObjSuf        = o
SrcSuf        = cxx
HeadSuf       = h
DllSuf        = so
DictpSuf      = Dict
DictExt       = h
DictSuf       = $(DictpSuf).$(DictExt)
DictObjSuf    = $(DictpSuf).$(ObjSuf)
OutPutOpt     = -o

#CXX           = g++
CXXFLAGS      = -O2 -fPIC -g -Wall
LIBS          = 
LDFLAGS       = -O
SOFLAGS       = -shared

CXX = g++
CC = gcc
LD = g++

ifeq ($(HYPHI_ENV),3.4.6)
CXX = g++-3.4
CC = gcc-3.4
LD = g++-3.4
#$(shell echo "---> gcc 3.4.4 will be used"):
else
CXX = g++
CC = gcc
LD = g++
endif

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
MACOSXTARGET := MACOSX_DEPLOYMENT_TARGET=10.$(MACOSX_MINOR)
ifeq ($(MACOSX_MINOR),5)
MACOSX_MINOR  = 4
endif
CXX           = c++
CXXFLAGS      = -O0 -fPIC -g
#-pipe -Wall -W -Woverloaded-virtual
LD            = $(MACOSXTARGET) c++
LDFLAGS       = -O -bind_at_load
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
UNDEFOPT      = dynamic_lookup
ifneq ($(MACOSX_MINOR),4)
ifneq ($(MACOSX_MINOR),3)
UNDEFOPT      = suppress
LD            = c++
endif
endif
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)
endif


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS) $(ROOTGLIBS) -lPhysics 

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : all clean clean_all

DEPEND        = makedepend $(CXXFLAGS)

#--------------------------------------------------------------------
Headers = TG4Sol_Event.hh TG4Sol_Hit.hh

Sources = TG4Sol_Event.cc TG4Sol_Hit.cc 

DictObj = DictionaryEventG4Sol.$(DictObjSuf)

ObjectEvent = TG4Sol_Event.$(ObjSuf) TG4Sol_Hit.$(ObjSuf)

G4HypEventSo = libEventG4Sol.$(DllSuf)

#-------------------------------------------

all: libG4_Event

libG4_Event : $(G4HypEventSo)

$(G4HypEventSo) : $(ObjectEvent) $(DictObj)
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
#ln -sf $(FromBindir)/$(G4HypEventSo) $(Bindir)/$(G4HypEventSo)
	@echo "---> $@ done"

clean:
	rm -f *.$(DictpSuf).* *.o *.so *.dylib

clean_all : clean

print :
	@echo "################################################################################"
	@echo "#                                                                               "
	@echo "# --> set up from HypHiLogin : $(HYPHI_ENV) / compiler selection : $(CXX) $(CC) "
	@echo "#                                                                               "
	@echo "################################################################################"

#--------------------------------------------------------------------
.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

$(DictObj):
	rm -f $*.*; \
	    rootcint $*.cc -s $(G4HypEventSo) -c $(Headers) TG4Sol_LinkDef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc

#.$(HeadSuf).$(DictObjSuf):
#	rm -f $*.$(DictpSuf).*; \
#	    rootcint $*.$(DictpSuf).cc -c $<; \
#	    $(CXX) $(CXXFLAGS) -c $*.$(DictpSuf).cc
