Bindir = ..
FromBindir = src
FromLibdir = ../src
DEPDIR = .deps
df = $(DEPDIR)/$(*F)


ifndef HYPHI_ENV
HYPHI_ENV := $(shell g++ --version | awk '/[0-9].[0-9].[0-9]/ { print $3 }' )
#  echo "HYPHI_ENV set up to $(HYPHI_ENV)"
endif

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)
ROOTVERSION  := $(shell root-config --version | sed -n "s/\([0-9]\)\.[0-9]*\/[0-9]*/\1/p")

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
CXXFLAGS      = -O3 -fPIC -g -std=c++1y #-Wall -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
LIBS          = -L./Genfit2_v2/lib -L./Ana_Event -L./EventG4Sol -L./FairBase -L.
LDFLAGS       = -O
SOFLAGS       = -shared

CXX = g++
CC =  gcc
LD = g++


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

GENFITV2FLAGS = -I./Genfit2_v2/include


CXXFLAGS     += $(ROOTCFLAGS) $(GENFITV2FLAGS) -I/opt/local/include #-I./boost
LIBS         += $(ROOTLIBS) $(ROOTGLIBS) -lEG -lGeom -lPhysics -lMinuit -lEventG4Sol -lTreePlayer -lHyphiFairBase -lMCAnaEventG4Sol

GENFIT_LIB = -lgenfit2
LIB_MERGER =  $(ROOTLIBS) -lAna_Event -lAna_Event_Out
LIB_PROFILING = #-L/home/christophe/Workspace/gperftools-2.4/lib -ltcmalloc #-lprofiler

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : clean libGenfit2_v2.so libEventG4Sol.so print clean_all check-syntax

DEPEND        = makedepend $(CXXFLAGS)

#--------------------------------------------------------------------
Sources = MainSimu.cc FullRecoEvent.cc Ana_Hist.cc  
Headers = $(patsubst %.cc,%.hh,$(Sources))
Objects = $(patsubst %.cc,%.o,$(Sources))

Sources_Proc =  TBuildDetectorLayerPlaneDAF.cc CheckField.cc #TBayesFinder.cc 

Headers_Proc = $(patsubst %.cc,%.h,$(Sources_Proc))
Objects_Proc = $(patsubst %.cc,%.o,$(Sources_Proc))

Sources_Proc_genfit_dep = THyphiAttributes.cc FullRecoTask.cc TKalmanFilter_DAF.cc 
Objects_Proc_genfit_dep = $(patsubst %.cc,%.o,$(Sources_Proc_genfit_dep))

CHK_SOURCES += $(Sources) $(Sources_Proc) $(Sources_Proc_genfit_dep) 
CHK_HEADERS = $(patsubst %.cc,%.hh,$(CHK_SOURCES))

#Objects_Proc_old_genfit = $(patsubst %.cc,%.o,$(Sources_Proc_old_genfit))

# Mathsources = MathematicalTools.cc
# MathObj = MathematicalTools.o
# MathSo = libHyphiMathTools.$(DllSuf)

all: print $(MathSo) libEventG4Sol.so libFairBase.so libGenfit2_v2.so libMCEvent.so $(Bindir)/MainSimu

$(Bindir)/MainSimu: $(Objects) $(Objects_Proc) $(Objects_Proc_genfit_dep) $(MathSo) libEventG4Sol.so libFairBase.so libGenfit2_v2.so
	$(CXX) $(CXXFLAGS)  -o $(Bindir)/MainSimu $(Objects) $(Objects_Proc) $(Objects_Proc_genfit_dep) $(LIBS) $(GENFIT_LIB) $(LIB_PROFILING)

# $(Bindir)/BatchMerger: BatchMerger.o libAnaEventOut.so libAna_Event.so
# 	$(CXX) $(CXXFLAGS) BatchMerger.o $(LIBS) $(LIB_MERGER) -o $(Bindir)/BatchMerger



libGenfit2_v1.so :
	make -C Genfit2_v1/build
	make -C Genfit2_v1/build install
	ln -sf $(FromLibdir)/Genfit2_v1/lib/libgenfit2.$(DllSuf).2.0 $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/Genfit2_v1/lib/coreRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v1/lib/eventDisplayRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v1/lib/fittersRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v1/lib/GBLRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v1/lib/measurementsRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v1/lib/trackRepsRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v1/lib/utilitiesRootDict_rdict.pcm $(Bindir)/lib/
endif
	@echo "---> $@ done"

libGenfit2_v2.so :
	make -C Genfit2_v2/build
	make -C Genfit2_v2/build install
	ln -sf $(FromLibdir)/Genfit2_v2/lib/libgenfit2.$(DllSuf).2.1 $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/Genfit2_v2/lib/coreRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v2/lib/eventDisplayRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v2/lib/fittersRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v2/lib/GBLRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v2/lib/measurementsRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v2/lib/trackRepsRootDict_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Genfit2_v2/lib/utilitiesRootDict_rdict.pcm $(Bindir)/lib/
endif
	@echo "---> $@ done"



libFairBase.so :
	make -C FairBase
	ln -sf $(FromLibdir)/FairBase/libHyphiFairBase.$(DllSuf) $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/FairBase/libHyphiFairBase_rdict.pcm $(Bindir)/lib/
endif
ifeq ($(ARCH),macosx)
	ln -sf $(FromLibdir)/FairBase/libHyphiFairBase.$(MacDllSuf) $(Bindir)/lib/
endif
	@echo "---> $@ done"


libEventG4Sol.so :
	make -C EventG4Sol
	ln -sf $(FromLibdir)/EventG4Sol/libEventG4Sol.$(DllSuf) $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/EventG4Sol/libEventG4Sol_rdict.pcm $(Bindir)/lib/
endif
ifeq ($(ARCH),macosx)
	ln -sf $(FromLibdir)/EventG4Sol/libEventG4Sol.$(MacDllSuf) $(Bindir)/lib/
endif
	@echo "---> $@ done"

libMCEvent.so :
	make -C Ana_Event
	ln -sf $(FromLibdir)/Ana_Event/libMCAnaEventG4Sol.$(DllSuf) $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/Ana_Event/libMCAnaEventG4Sol_rdict.pcm $(Bindir)/lib/
endif
ifeq ($(ARCH),macosx)
	ln -sf $(FromLibdir)/Ana_Event/libMCAnaEvent.$(MacDllSuf) $(Bindir)/lib/
endif
	@echo "---> $@ done"


# $(MathSo) : $(MathObj)
# ifeq ($(PLATFORM),macosx)
# # We need to make both the .dylib and the .so
# 	$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
# ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
# ifeq ($(MACOSX_MINOR),4)
# 	ln -sf $@ $(subst .$(DllSuf),.so,$@)
# else
# 	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
# endif
# endif
# else
# 	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
# endif
# 	ln -sf $(FromLibdir)/$@ $(Bindir)/lib/
# 	@echo "---> $@ done"


clean:
	rm -f $(Bindir)/BatchMerger $(Bindir)/MainSimu $(Bindir)/*.so $(Bindir)/lib/*.so $(Bindir)/lib/*.pcm $(Bindir)/lib/*.dylib *.$(DictpSuf).* *.o *.so *.dylib *pcm

clean_all: clean
	make -C Genfit2_v2/build clean
	make -C EventG4Sol clean_all
	make -C FairBase clean_all
	make -C Ana_Event clean_all
print :
	@echo "################################################################################"
	@echo "#                                                                               "
	@echo "# --> set up from HypHiLogin : $(HYPHI_ENV) / compiler selection : $(CXX) $(CC) "
	@echo "#                                                                               "
	@echo "################################################################################"

check-syntax:
#	$(CXX) $(CXXFLAGS) -Wall -Wextra -fsyntax-only $(CHK_SOURCES)
	clang-tidy-4.0 $(CHK_SOURCES) $(CHK_HEADERS) -config='' -- -std=c++14 -I/usr/include/c++/5/ -I/home/christophe/root-6/include -I./Genfit2_v2/include
#depend:
#	 makedepend -f- -- $(CXXFLAGS) -- *.cc > Fullrec.depends
#-include Fullrec.depends

#--------------------------------------------------------------------
#.$(SrcSuf).$(ObjSuf):
#%.$(ObjSuf) : %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -c $<

#.$(HeadSuf).$(DictObjSuf):
#	rm -f $*.$(DictpSuf).*; \
#	    rootcint $*.$(DictpSuf).cc -c $<; \
#	    $(CXX) $(CXXFLAGS) -c $*.$(DictpSuf).cc

#%.$(DictObjSuf): %.$(HeadSuf)
#	rm -f $*.$(DictpSuf).*; \
#	    rootcint $*.$(DictpSuf).cc -c $<; \
#	    $(CXX) $(CXXFLAGS) -c $*.$(DictpSuf).cc

#$(Objects_Proc_genfit_dep): %.$(ObjSuf): %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) $(GENFITFLAGS) -c $<

%.$(ObjSuf) : %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -MD -c $<
	@cp $*.d $(df).P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
            rm -f $*.d

-include $(Sources:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(Sources_Proc:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(Sources_Proc_genfit_dep:%.$(SrcSuf)=$(DEPDIR)/%.P)
#-include $(Sources_Proc_old_genfit:%.$(SrcSuf)=$(DEPDIR)/%.P)


%.$(DictObjSuf): %.$(HeadSuf)
	rm -f $*.$(DictpSuf).*; \
	    rootcint $*.$(DictpSuf).cc -c $<; \
	    $(CXX) $(CXXFLAGS) -c $*.$(DictpSuf).cc

$(Objects_Proc_genfit_dep): %.$(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) $(GENFITV2FLAGS) -MD -c $<
	@cp $*.d $(df).P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
            rm -f $*.d

$(Objects_Second): %.$(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) $(GENFITFLAGS) -MD -c $<
	@cp $*.d $(df).P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
            rm -f $*.d

$(Objects_Simulator): %.$(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) $(GENFITFLAGS) -MD -c $<
	@cp $*.d $(df).P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
            rm -f $*.d

$(Dic_AliKF):
	rm -f $*.*; \
	    rootcint $*.cc -c $(Headers_AliKF) AliKF_linkdef.h; \
	    $(CXX) $(CXXFLAGS) -c $*.cc
