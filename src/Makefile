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
CXXFLAGS      = -Wall -O3 -march=native -fPIC -g -DSPDLOG_COMPILED_LIB -DZMQ_BUILD_DRAFT_API -DNonhomogeneousField # -fsanitize=address -fno-omit-frame-pointer #-Wall -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
#CXXFLAGS     += -fsanitize=memory -fno-omit-frame-pointer
#CXXFLAGS     += -fsanitize=undefined -fno-omit-frame-pointer
#CXXFLAGS      += -fsanitize=leak -fsanitize-blacklist=$(shell root-config --etcdir)/lsan-root.supp
LIBS          = -L./TrickTrack/lib -Lmsgpack-c -Llibzmq/lib -Lspdlog/lib/spdlog -L./Genfit2_v2/lib -L./Ana_Event -L./EventG4Sol -L./EventWASAUnpack -L./HitAna -L./KFParticle/install/lib -L./FairBase -L.

LDFLAGS       = '-Wl,-rpath,$$ORIGIN/lib'
SOFLAGS       = -shared

CXX = clang++
CC =  clang
LD = clang++


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

GEANTCFLAGS := $(shell geant4-config --cflags)
GEANTLIBS   := $(shell geant4-config --libs)

GENFITV2FLAGS = -I./Genfit2_v2/include
SPDLOGFLAGS = -I./spdlog/include
ZMQFLAGS = -I./cppzmq/include -I./libzmq/include
MSGFLAGS = -I./msgpack-c/include
TRICKTRACKFLAGS = -I./TrickTrack/include
KFPARTICLEFLAGS = -I./KFParticle/install/include
ORMFLAGS = -I./sqlite_orm/include

CXXFLAGS     += $(MSGFLAGS) $(ZMQFLAGS) $(SPDLOGFLAGS) $(ROOTCFLAGS) $(GEANTCFLAGS) $(GENFITV2FLAGS) $(TRICKTRACKFLAGS) $(KFPARTICLEFLAGS) $(ORMFLAGS) -I/opt/local/include -I/usr/include/eigen3 #-I./boost
LIBS_SIM         = $(LIBS) $(ROOTLIBS) $(ROOTGLIBS) -lEG -lGeom -lPhysics -lMinuit -lEventG4Sol -lEventWASAUnpack -lHitAna -lTreePlayer -lHyphiFairBase -lMCAnaEventG4Sol -lAna_Event -lsqlite3 $(GEANTLIBS)
LIBS_DATA        = $(LIBS) $(ROOTLIBS) $(ROOTGLIBS) -lEG -lGeom -lPhysics -lMinuit -lEventWASAUnpack -lHitAna -lEventG4Sol -lTreePlayer -lHyphiFairBase -lAna_Event -lMCAnaEventG4Sol -lsqlite3 $(GEANTLIBS)

GENFIT_LIB = -lgenfit2
SPDLOG_LIB = -lspdlog
ZMQ_LIB = -lzmq
MSG_LIB = -lmsgpackc
TT_LIB = -ltricktrack
KF_LIB = -lKFParticle

LIB_MERGER =  $(ROOTLIBS) -lAna_Event -lAna_Event_Out
LIB_PROFILING = #-L/home/christophe/Workspace/gperftools-2.4/lib -ltcmalloc #-lprofiler

.SUFFIXES : .$(HeadSuf) .$(DictSuf) .$(ObjSuf) .$(SrcSuf) .$(DictObjSuf)
.PHONY : clean libGenfit2_v2.so libEventG4Sol.so print clean_all check-syntax

DEPEND        = makedepend $(CXXFLAGS)

#--------------------------------------------------------------------
SourcesSimu = MainSimu.cc
Objects_Simu = $(patsubst %.cc,%.o,$(SourcesSimu))

SourcesData = MainData.cc  
Objects_Data = $(patsubst %.cc,%.o,$(SourcesData))

Sources_Comm = FullRecoEvent.cc FullRecoEventZMQ.cc FullRecoConfig.cc Ana_Hist.cc  
Headers_Comm = $(patsubst %.cc,%.hh,$(Sources_Comm))
Objects_Comm = $(patsubst %.cc,%.o,$(Sources_Comm))
 
Sources_Proc =  TBuildDetectorLayerPlaneDAF.cc TBuildWASACalibrationLayerPlane.cc TBuildRestarter.cc CheckField.cc TBuildDetectorLayerPlaneDAF_MT.cc TMergerOutput_MT.cc TBuildDetectorLayerPlaneDAF_ZMQ.cc TMergerOutput_ZMQ.cc TBayesFinder.cc TCheckRZ.cc TFlatMCOutputML.cc TPrimaryVertex.cc TDecayVertex.cc TRiemannFinder.cc TFindingPerf.cc TCheckFiberTrack.cc #TFinderCM.cc

Headers_Proc = $(patsubst %.cc,%.h,$(Sources_Proc))
Objects_Proc = $(patsubst %.cc,%.o,$(Sources_Proc))

Sources_Proc_genfit_dep = THyphiAttributes.cc FullRecoTask.cc TKalmanFilter_DAF.cc FullRecoTaskMT.cc TKalmanFilter_DAF_MT.cc
Objects_Proc_genfit_dep = $(patsubst %.cc,%.o,$(Sources_Proc_genfit_dep))

Source_Proc_genfit_ZMQ_dep = FullRecoTaskZMQ.cc TKalmanFilter_DAF_ZMQ.cc
Objects_Proc_genfit_ZMQ_dep = $(patsubst %.cc,%.o,$(Source_Proc_genfit_ZMQ_dep))

CHK_SOURCES += $(Sources_Comm) $(Sources_Proc) $(Sources_Proc_genfit_dep) $(Sources_Proc_genfit_ZMQ_dep)
CHK_HEADERS = $(patsubst %.cc,%.hh,$(CHK_SOURCES))

#Objects_Proc_old_genfit = $(patsubst %.cc,%.o,$(Sources_Proc_old_genfit))

# Mathsources = MathematicalTools.cc
# MathObj = MathematicalTools.o
# MathSo = libHyphiMathTools.$(DllSuf)

all: print $(MathSo) libEventG4Sol.so libEventWASAUnpack.so libHitAna.so libFairBase.so libGenfit2_v2.so libspdlog.a libMCEvent.so libzmq.so libmsgpack.so $(Bindir)/MainSimu $(Bindir)/MainData $(Bindir)/MainMonitor 

$(Bindir)/MainSimu: $(Objects_Simu) $(Objects_Comm) $(Objects_Proc) $(Objects_Proc_genfit_dep) $(Objects_Proc_genfit_ZMQ_dep) $(MathSo) libEventG4Sol.so libFairBase.so libGenfit2_v2.so libspdlog.a libtricktrack.so libKFParticle.so
	$(CXX) $(CXXFLAGS) -o $(Bindir)/MainSimu $(Objects_Simu) $(Objects_Comm) $(Objects_Proc) $(Objects_Proc_genfit_dep) $(Objects_Proc_genfit_ZMQ_dep) $(LIBS_SIM) $(GENFIT_LIB) $(LIB_PROFILING) $(SPDLOG_LIB) $(ZMQ_LIB) $(MSG_LIB) $(TT_LIB) $(KF_LIB) $(LDFLAGS)

$(Bindir)/MainData: $(Objects_Data) $(Objects_Comm) $(Objects_Proc) $(Objects_Proc_genfit_dep) $(Objects_Proc_genfit_ZMQ_dep) $(MathSo) libEventWASAUnpack.so libHitAna.so libFairBase.so libGenfit2_v2.so libspdlog.a libtricktrack.so libKFParticle.so
	$(CXX) $(CXXFLAGS) -o $(Bindir)/MainData $(Objects_Data) $(Objects_Comm) $(Objects_Proc) $(Objects_Proc_genfit_dep) $(Objects_Proc_genfit_ZMQ_dep) $(LIBS_DATA) $(GENFIT_LIB) $(LIB_PROFILING) $(SPDLOG_LIB) $(ZMQ_LIB) $(MSG_LIB) $(TT_LIB) $(KF_LIB) $(LDFLAGS) 

$(Bindir)/MainMonitor: MainMonitor.o FullRecoConfig.o
	$(CXX) $(CXXFLAGS) -o $(Bindir)/MainMonitor MainMonitor.o FullRecoConfig.o $(LIBS_SIM) -lRHTTP $(LIB_PROFILING) $(SPDLOG_LIB) $(ZMQ_LIB) $(MSG_LIB) $(LDFLAGS)

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
	ln -sf $(FromLibdir)/Genfit2_v2/lib/libgenfit2.$(DllSuf).2.2 $(Bindir)/lib/
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

libspdlog.a :
	make -C spdlog/build
	make -C spdlog/build install

libtricktrack.so :
	make -C TrickTrack/build
	make -C TrickTrack/build install
	ln -sf $(FromLibdir)/TrickTrack/lib/libtricktrack.$(DllSuf) $(Bindir)/lib/
	@echo "---> $@ done"

libKFParticle.so :
	make -C KFParticle/build
	make -C KFParticle/build install
	ln -sf $(FromLibdir)/KFParticle/install/lib/libKFParticle.$(DllSuf) $(Bindir)/lib/
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

libEventWASAUnpack.so :
	make -C EventWASAUnpack
	ln -sf $(FromLibdir)/EventWASAUnpack/libEventWASAUnpack.$(DllSuf) $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/EventWASAUnpack/libEventWASAUnpack_rdict.pcm $(Bindir)/lib/
endif
ifeq ($(ARCH),macosx)
	ln -sf $(FromLibdir)/EventWASAUnpack/libEventWASAUnpack.$(MacDllSuf) $(Bindir)/lib/
endif
	@echo "---> $@ done"

libHitAna.so :
	make -C HitAna
	ln -sf $(FromLibdir)/HitAna/libHitAna.$(DllSuf) $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/HitAna/libHitAna_rdict.pcm $(Bindir)/lib/
endif
ifeq ($(ARCH),macosx)
	ln -sf $(FromLibdir)/HitAna/libHitAna.$(MacDllSuf) $(Bindir)/lib/
endif
	@echo "---> $@ done"

libMCEvent.so :
	make -C Ana_Event
	ln -sf $(FromLibdir)/Ana_Event/libMCAnaEventG4Sol.$(DllSuf) $(Bindir)/lib/
	ln -sf $(FromLibdir)/Ana_Event/libAna_Event.$(DllSuf) $(Bindir)/lib/
	ln -sf $(FromLibdir)/Ana_Event/libAna_WasaEvent.$(DllSuf) $(Bindir)/lib/
ifeq ($(ROOTVERSION),6)
	ln -sf $(FromLibdir)/Ana_Event/libMCAnaEventG4Sol_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Ana_Event/libAna_WasaEvent_rdict.pcm $(Bindir)/lib/
	ln -sf $(FromLibdir)/Ana_Event/libAna_WasaEvent_rdict.pcm $(Bindir)/lib/libAna_Event_rdict.pcm
endif
ifeq ($(ARCH),macosx)
	ln -sf $(FromLibdir)/Ana_Event/libMCAnaEvent.$(MacDllSuf) $(Bindir)/lib/
endif
	@echo "---> $@ done"

libzmq.so :
	ln -sf $(FromLibdir)/libzmq/lib/libzmq.so.5 $(Bindir)/lib/
	@echo "---> $@ done"

libmsgpack.so :
	ln -sf $(FromLibdir)/msgpack-c/libmsgpackc.so.2 $(Bindir)/lib/
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
	rm -f $(Bindir)/BatchMerger $(Bindir)/MainMonitor $(Bindir)/MainSimu $(Bindir)/*.so $(Bindir)/lib/*.so $(Bindir)/lib/*.pcm $(Bindir)/lib/*.dylib *.$(DictpSuf).* *.o *.so *.dylib *pcm

clean_all: clean
	make -C Genfit2_v2/build clean
	make -C EventG4Sol clean_all
	make -C FairBase clean_all
	make -C Ana_Event clean_all
	make -C EventWASAUnpack clean_all
	make -C HitAna clean_all
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

-include $(SourcesSimu:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(SourcesData:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(Sources_Comm:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(Sources_Proc:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(Sources_Proc_genfit_dep:%.$(SrcSuf)=$(DEPDIR)/%.P)
-include $(Sources_Proc_genfit_ZMQ_dep:%.$(SrcSuf)=$(DEPDIR)/%.P)
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

$(Objects_Proc_genfit_ZMQ_dep): %.$(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) $(GENFITV2FLAGS) $(ZMQFLAGS) $(MSGFLAGS) -MD -c $<
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
