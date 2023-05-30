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
CXXFLAGS      = -O2 -march=native -fPIC -g -Wall
LIBS          = 
LDFLAGS       = -O
SOFLAGS       = -shared

CXX = g++
CC = gcc
LD = g++


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
Headers = WASAUnpackBranch.hh

Sources = WASAUnpackBranch.cxx

DictObj = DictionaryEventWASAUnpack.$(DictObjSuf)

ObjectEvent = WASAUnpackBranch.$(ObjSuf)

EventWASAUnpackSo = libEventWASAUnpack.$(DllSuf)

#-------------------------------------------

all: libWASAUnpack_Event

libWASAUnpack_Event : $(EventWASAUnpackSo)

$(EventWASAUnpackSo) : $(ObjectEvent) $(DictObj)
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
#ln -sf $(FromBindir)/$(EventWASAUnpackSo) $(Bindir)/$(EventWASAUnpackSo)
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
	    rootcint $*.cxx -s $(EventWASAUnpackSo) $(Headers) WASAUnpack_LinkDef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cxx

#.$(HeadSuf).$(DictObjSuf):
#	rm -f $*.$(DictpSuf).*; \
#	    rootcint $*.$(DictpSuf).cc -c $<; \
#	    $(CXX) $(CXXFLAGS) -c $*.$(DictpSuf).cc
