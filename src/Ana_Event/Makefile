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

ROOTDICGEN = rootcint
ifeq ($(ROOTVERSION),6)
ROOTDICGEN = rootcling
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
ROOTLIBSDIR  := $(shell root-config --libdir)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS         +=-L$(ROOTLIBSDIR) -lPhysics 
#$(ROOTLIBS) $(ROOTGLIBS) -lPhysics 

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : all clean clean_all

DEPEND        = makedepend $(CXXFLAGS)

#--------------------------------------------------------------------

EventH_v14 = Ana_EventNew_v14.hh
EventH_Com_v14 = THypernucleus.hh THyphiHitDet.hh THyphiTrack_v2.hh TMcParticle.hh
EventH_v16 = Ana_EventNew_v16.hh
EventH_Com_v16 = THypernucleus.hh THyphiHitDet.hh THyphiTrack_v3.hh TMcParticle.hh
MCEventH = MCAnaEvent.hh
MCEventH_Com = TMcHit.hh TMcParticle.hh THyphiTrack_v3.hh

MCEventG4SolH = MCAnaEventG4Sol.hh
MCEventG4SolH_Com = TMcHit.hh TMcParticle.hh THyphiTrack_v4.hh


EventSrc_v14 = $(patsubst %.hh,%.cc,$(EventH_v14))
EventSrc_Com_v14 = $(patsubst %.hh,%.cc,$(EventH_Com_v14))
EventSrc_v16 = $(patsubst %.hh,%.cc,$(EventH_v16))
EventSrc_Com_v16 = $(patsubst %.hh,%.cc,$(EventH_Com_v16))

MCEventSrc = $(patsubst %.hh,%.cc,$(MCEventH))
MCEventSrc_Com = $(patsubst %.hh,%.cc,$(MCEventH_Com))

MCEventG4SolSrc = $(patsubst %.hh,%.cc,$(MCEventG4SolH))
MCEventG4SolSrc_Com = $(patsubst %.hh,%.cc,$(MCEventG4SolH_Com))


EventDict_Com = THypernucleus.$(DictObjSuf) THyphiSimpleHit.$(DictObjSuf) TTrackSimple.$(DictObjSuf) TMcParticle.$(DictObjSuf)
EventDict_v14 = Ana_EventNew_v14.$(DictObjSuf)
EventDict_v16 = Ana_EventNew_v16.$(DictObjSuf)

MCEventDict = MCAnaEvent.$(DictObjSuf)
MCEventG4SolDict = MCAnaEventG4Sol.$(DictObjSuf)

EventO_Com = THypernucleus.$(ObjSuf) THyphiSimpleHit.$(ObjSuf) TTrackSimple.$(ObjSuf) TMcParticle.$(ObjSuf)

EventO_v14 = $(patsubst %.hh,%.$(ObjSuf),$(EventH_v14))
EventO_Com_v14 = $(patsubst %.hh,%.$(ObjSuf),$(EventH_Com_v14))

EventO_v16 = $(patsubst %.hh,%.$(ObjSuf),$(EventH_v16))
EventO_Com_v16 = $(patsubst %.hh,%.$(ObjSuf),$(EventH_Com_v16))

MCEventO = $(patsubst %.hh,%.$(ObjSuf),$(MCEventH))
MCEventO_Com = $(patsubst %.hh,%.$(ObjSuf),$(MCEventH_Com))

MCEventG4SolO = $(patsubst %.hh,%.$(ObjSuf),$(MCEventG4SolH))
MCEventG4SolO_Com = $(patsubst %.hh,%.$(ObjSuf),$(MCEventG4SolH_Com))

EventSo_v14 = libAna_Event_v14.$(DllSuf)
EventSo_v16 = libAna_Event_v16.$(DllSuf)
MCEventSo = libMCAnaEvent.$(DllSuf)
MCEventG4SolSo = libMCAnaEventG4Sol.$(DllSuf)

EventSo_Commun_Name = libAna_Event

all: libAna_Event 

libAna_Event : $(MCEventG4SolSo) 
	ln -sf $(EventSo_v16) $(EventSo_Commun_Name).$(DllSuf)
ifeq ($(PLATFORM),macosx)
	ln -sf $(EventSo_v16) $(EventSo_Commun_Name).$(MacDllSuf)
endif


$(EventSo_v14): $(EventO_v14) $(EventO_Com_v14) $(EventDict_v14) $(EventDict_Com_v14)
# $(EventSo_v12): $(EventO_v12) $(EventO_Com_v12) $(EventDict_v12) $(EventDict_Com_v12)
#$(EventSo_v10): $(EventO_v10) $(EventO_Com_v2) $(EventDict_v10) $(EventDict_Com_v2)
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
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(LIBS)
endif
#	ln -sf $(FromLibdir)/$(EventSo_v10) $(Bindir)/lib
	@echo "---> $@ done"

$(EventSo_v16): $(EventO_v16) $(EventO_Com_v16) $(EventDict_v16) $(EventDict_Com_v16)
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
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(LIBS)
endif
#	ln -sf $(FromLibdir)/$(EventSo_v10) $(Bindir)/lib
	@echo "---> $@ done"



$(MCEventSo): $(MCEventO) $(MCEventO_Com) $(MCEventDict)
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
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(LIBS)
endif
#	ln -sf $(FromLibdir)/$(EventSo_v10) $(Bindir)/lib
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

 .$(HeadSuf).$(DictObjSuf):
	rm -f $*.$(DictpSuf).*; \
	    $(ROOTDICGEN) $*.$(DictpSuf).cc $<; \
	    $(CXX) $(CXXFLAGS) -c $*.$(DictpSuf).cc

$(EventDict_v14):
	rm -f $*.*; \
	    $(ROOTDICGEN) $*.cc -s $(EventSo_v14) $(EventH_v14) $(EventH_Com_v14) Ana_EventNew_v14_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc

$(EventDict_v16):
	rm -f $*.*; \
	    $(ROOTDICGEN) $*.cc -s $(EventSo_v16) $(EventH_v16) $(EventH_Com_v16) Ana_EventNew_v16_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc

$(MCEventDict):
	rm -f $*.*; \
	    $(ROOTDICGEN) $*.cc -s $(MCEventSo) $(MCEventH) $(MCEventH_Com) MCAnaEvent_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc

$(MCEventG4SolDict):
	rm -f $*.*; \
	    $(ROOTDICGEN) $*.cc -s $(MCEventG4SolSo) $(MCEventG4SolH) $(MCEventG4SolH_Com) MCAnaEventG4Sol_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc
