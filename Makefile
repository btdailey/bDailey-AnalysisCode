#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES     = -I/rh5stuff/64bit/include -I/usr/include
SYSLIBS         = -L/rh5stuff/64bit/lib -L/usr/lib64

ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR)
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
else
ANITA_UTIL_LIB_DIR=/usr/local/lib
ANITA_UTIL_INC_DIR=/usr/local/include
ifdef EVENT_READER_DIR
LD_ANITA_UTIL=-L$(EVENT_READER_DIR)
INC_ANITA_UTIL=-I$(EVENT_READER_DIR)
endif
endif

#Toggles the FFT functions on and off
USE_FFT_TOOLS=1

ifdef USE_FFT_TOOLS
#FFTLIBS = -lRootFftwWrapper -lfftw3
FFTLIBS = -lRootFftwWrapper -lfftw3
FFTFLAG = -DUSE_FFT_TOOLS
else
FFTLIBS =
FFTFLAG =
endif

#Generic and Site Specific Flags
CXXFLAGS     += -g -pg $(ROOTCFLAGS) $(FFTFLAG) $(INC_ANITA_UTIL) $(SYSINCLUDES) -I$(ANITA_ANALYSIS_HOOVER) -I$(EVENT_SIMULATION_DIR) -I/home/dailey.110/Healpix_3.30/include -I/data/anita/btdailey/Healpix_stuff/cfitsio/include -I/home/dailey.110/Healpix_3.30/src/cxx/generic_gcc/include
LDFLAGS      += -g $(ROOTLDFLAGS) 

LIBS          = -L/home/dailey.110/analysis/ -L/data/anita/btdailey/Healpix_stuff/cfitsio/lib -L/home/dailey.110/Healpix_3.30/lib -L/home/dailey.110/Healpix_3.30/src/cxx/generic_gcc/lib $(ROOTLIBS) -lchealpix -lcfitsio -lMathMore -lMinuit -L/data/anita/btdailey/analysis_info/data $(SYSLIBS) $(LD_ANITA_UTIL) $(FFTLIBS) -L$(ANITA_ANALYSIS_HOOVER)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling
ROOT_LIBRARY = libMyCorrelator.${DLLSUF}
LIB_OBJS =MyCorrelator.o corrDict.o
CLASS_HEADERS = MyCorrelator.h
#EXE = ssEventMaker

all : $(ROOT_LIBRARY)


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DLLSUF),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	 $(OutPutOpt) $(subst .$(DLLSUF),.so,$@)
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) $(LIBS) $(GLIBS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


corrDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(INC_ANITA_UTIL) -I$(ANITA_ANALYSIS_HOOVER) -I$(EVENT_SIMULATION_DIR) $(CLASS_HEADERS) LinkDef.h

install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	cp $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
else
	cp $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	cp  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)

#executable: $(EXE)

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
	@rm -f $(EXE)
#############################################################################



