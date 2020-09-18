Bindir    = .

ifndef HYPHI_ENV
HYPHI_ENV := $(shell g++ --version | awk '/[0-9].[0-9].[0-9]/ { print $3 }' )
#  echo "HYPHI_ENV set up to $(HYPHI_ENV)"
endif

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

ifeq ($(ARCH),macosx64)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
MACOSXTARGET := MACOSX_DEPLOYMENT_TARGET=10.$(MACOSX_MINOR)
ifeq ($(MACOSX_MINOR),6)
MACOSX_MINOR  = 4
endif
CXX           = g++
CXXFLAGS      = -O2 -fPIC -g
#-pipe -Wall -W -Woverloaded-virtual
LD            = $(MACOSXTARGET) g++
LDFLAGS       = -O -bind_at_load
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
UNDEFOPT      = dynamic_lookup
ifneq ($(MACOSX_MINOR),4)
ifneq ($(MACOSX_MINOR),3)
UNDEFOPT      = suppress
LD            = g++
endif
endif
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)
endif


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS) $(ROOTGLIBS) -lGeom 
LIBS_EVE     += -lEve

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : all clean clean_all

DEPEND        = makedepend $(CXXFLAGS)




Sources = FairField.cxx HypFieldMap.cxx HypFieldMapData.cxx HypFieldMapFull.cxx HypFieldMapFullData.cxx HypConstField.cxx FrsSksHypFieldMapFull.cxx FrsSksHypFieldMapFullData.cxx FrsSolenoidHypField.cxx
SourcesEve = EveAladinField.cxx EveSksField.cxx

Headers = $(patsubst %.cxx,%.h,$(Sources))
HeadersEve = $(patsubst %.cxx,%.h,$(SourcesEve))
Objects = $(patsubst %.cxx,%.o,$(Sources))
ObjectsEve = $(patsubst %.cxx,%.o,$(SourcesEve))


DictObj = FairBaseDictionary.$(DictObjSuf)
DictEveObj = EveFairBaseDictionary.$(DictObjSuf)

FairBase = libHyphiFairBase.$(DllSuf)
EveBase = libHyphiEve.$(DllSuf)

all: libHyphiFairBase libHyphiEve

libHyphiFairBase : $(FairBase)
libHyphiEve : $(EveBase)

$(FairBase): $(Objects) $(DictObj)
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
#ln -sf $(FromBindir)/$(Event_Phys) $(Bindir)/$(Event_Phys)
	@echo "---> $@ done"


$(EveBase): $(Objects) $(ObjectsEve) $(DictEveObj)
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
	@echo "---> $@ done"



clean:
	rm -f $(Bindir)/$(FairBase) *.$(DictpSuf).* *.o *.so *.pcm

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
	    rootcint $*.cc -s $(FairBase) $(Headers) FairBaseLinkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc

$(DictEveObj):
	rm -f $*.*; \
	    rootcint $*.cc -s $(EveBase) $(Headers) $(HeadersEve) EveFairLinkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc
