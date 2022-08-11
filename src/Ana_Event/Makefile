Bindir = .

ifndef HYPHI_ENV
HYPHI_ENV := $(shell g++ --version | awk '/[0-9].[0-9].[0-9]/ { print $3 }' )
#  echo "HYPHI_ENV set up to $(HYPHI_ENV)"
endif

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)
ROOTVERSION := $(shell root-config --version | cut -d . -f 1)

ObjSuf        = o
SrcSuf        = cc
HeadSuf       = hh
DllSuf        = so
MacDllSuf     = so
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

ROOTDICGEN = rootcint
ifeq ($(ROOTVERSION),6)
ROOTDICGEN = rootcling
endif


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTLIBSDIR  := $(shell root-config --libdir)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS         +=-L$(ROOTLIBSDIR) -lPhysics 

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : all clean clean_all

DEPEND        = makedepend $(CXXFLAGS)

#--------------------------------------------------------------------

EventH_Wasa = Ana_WasaEvent.hh AnaEvent_Metadata.hh
EventH_Com_Wasa = THyphiTrack_v4.hh TTrackCand.hh THypernucleus.hh TDataHit.hh

MCEventG4SolH = MCAnaEventG4Sol.hh AnaEvent_Metadata.hh
MCEventG4SolH_Com = TMcHit.hh TMcParticle.hh THyphiTrack_v4.hh TTrackCand.hh THypernucleus.hh


EventSrc_Wasa = $(patsubst %.hh,%.cc,$(EventH_Wasa))
EventSrc_Com_Wasa = $(patsubst %.hh,%.cc,$(EventH_Com_Wasa))

MCEventG4SolSrc = $(patsubst %.hh,%.cc,$(MCEventG4SolH))
MCEventG4SolSrc_Com = $(patsubst %.hh,%.cc,$(MCEventG4SolH_Com))

EventDict_Wasa = Ana_WasaEvent.$(DictObjSuf)

MCEventG4SolDict = MCAnaEventG4Sol.$(DictObjSuf)

EventO_Com = THypernucleus.$(ObjSuf) THyphiSimpleHit.$(ObjSuf) TTrackSimple.$(ObjSuf) TMcParticle.$(ObjSuf)

EventO_Wasa = $(patsubst %.hh,%.$(ObjSuf),$(EventH_Wasa))
EventO_Com_Wasa = $(patsubst %.hh,%.$(ObjSuf),$(EventH_Com_Wasa))

MCEventG4SolO = $(patsubst %.hh,%.$(ObjSuf),$(MCEventG4SolH))
MCEventG4SolO_Com = $(patsubst %.hh,%.$(ObjSuf),$(MCEventG4SolH_Com))

EventSo_Wasa = libAna_WasaEvent.$(DllSuf)
MCEventG4SolSo = libMCAnaEventG4Sol.$(DllSuf)

EventSo_Commun_Name = libAna_Event

all: libAna_Event $(MCEventG4SolSo)

libAna_Event : $(EventSo_Wasa)
	ln -sf $(EventSo_Wasa) $(EventSo_Commun_Name).$(DllSuf)

$(EventSo_Wasa): $(EventO_Wasa) $(EventO_Com_Wasa) $(EventDict_Wasa)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(LIBS)
	@echo "---> $@ done"


$(MCEventG4SolSo): $(MCEventG4SolO) $(MCEventG4SolO_Com) $(MCEventG4SolDict)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(LIBS)
	@echo "---> $@ done"


clean:
	rm -f *.$(DictpSuf).* *.o *.so *.dylib *.pcm

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


$(EventDict_Wasa):
	rm -f $*.*; \
	    $(ROOTDICGEN) $*.cc -s $(EventSo_Wasa) $(EventH_Wasa) $(EventH_Com_Wasa) Ana_WasaEvent_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc

$(MCEventG4SolDict):
	rm -f $*.*; \
	    $(ROOTDICGEN) $*.cc -s $(MCEventG4SolSo) $(MCEventG4SolH) $(MCEventG4SolH_Com) MCAnaEventG4Sol_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc
