Bindir = .

ifndef HYPHI_ENV
HYPHI_ENV := $(shell g++ --version | awk '/[0-9].[0-9].[0-9]/ { print $3 }' )
#  echo "HYPHI_ENV set up to $(HYPHI_ENV)"
endif

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)
ROOTVERSION := $(shell root-config --version | cut -d . -f 1)

ObjSuf        = o
SrcSuf        = cxx
HeadSuf       = hh
DllSuf        = so
MacDllSuf     = so
DictpSuf      = Dict
DictExt       = h
DictSuf       = $(DictpSuf).$(DictExt)
DictObjSuf    = $(DictpSuf).$(ObjSuf)
OutPutOpt     = -o

CXXFLAGS      = -O2 -march=native -fPIC -g -Wall
LIBS          = 
LDFLAGS       = -O
SOFLAGS       = -shared

CXX = g++
CC = gcc
LD = g++

ROOTDICGEN = rootcint
ifeq ($(ROOTVERSION),6)
ROOTDICGEN = rootcling
endif

WASAUNPACK = -I../EventWASAUnpack

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTLIBSDIR  := $(shell root-config --libdir)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS) $(WASAUNPACK)
LIBS         +=-L$(ROOTLIBSDIR) -lPhysics 

GENFITV2FLAGS = -I../Genfit2_v2/include
CXXFLAGS     += $(GENFITV2FLAGS)

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : all clean clean_all

DEPEND        = makedepend $(CXXFLAGS)

#--------------------------------------------------------------------

Headers = FiberAnalyzer.hh FiberHitAna.hh FiberHitXUV.hh FiberTrackAna.hh MDCHitAna.hh ParaManager.hh PSBEHitAna.hh PSBHitAna.hh PSFEHitAna.hh T0HitAna.hh TrackHit.hh MWDCHitAna.hh MWDCTracking.hh OpticsMom.hh S4SciHitAna.hh

Sources = $(patsubst %.$(HeadSuf),%.$(SrcSuf),$(Headers))

ObjectEvent = $(patsubst %.$(HeadSuf),%.$(ObjSuf),$(Headers))

HitAnaSo = libHitAna.$(DllSuf)

#-------------------------------------------

all: libAna_Hit

libAna_Hit : $(HitAnaSo)

$(HitAnaSo) : $(ObjectEvent)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(LIBS)
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

