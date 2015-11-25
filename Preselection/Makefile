BASEDIR = $(ROOTSORTERSYS)/core
TARGET   = main
MODULES  =  log theory analysis montecarlo data reconstruction detectors trackprocessing he3 analysiswrap
DICTMODULES = analysiswrap

CPPFLAGS := `$(BASEDIR)/bin/sorter-config -cpp` -g  
CXXFLAGS := `$(BASEDIR)/bin/sorter-config -cxx` -fno-strict-aliasing  -g -std=c++11
LDFLAGS  := `$(BASEDIR)/bin/sorter-config -ld -libs-wasa -libs-wasa-ana` -lWasaRecFD -lWasaRecSE -lWasaRecPS -lWasaRecMDC -lWasaRecCD -lWasaParameter -lWasaAnaRaw -lWasaRecFPC -lvt -g 

all: $(TARGET)

$(TARGET): $(TARGET).o 
	$(CXX) $(LDFLAGS) $^ -o $@

ifneq (xx,x$(MODULES)x)
$(TARGET): $(TARGET).so

$(TARGET).so: $(addsuffix .o,$(MODULES)) $(TARGET)Dict.o
	$(CXX) -fPIC -O2 -shared $^ -o $@

$(TARGET)Dict.cc: $(addsuffix .hh,$(DICTMODULES)) LinkDef.hh
	rootcint -f $@ -c $(CPPFLAGS) -p $^ 
endif

%.o: %.cc 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -p -c $< -o $@

%.so: %.cc
	rootcint -f $(subst .cc,Dict.cc,$<) -c $(CPPFLAGS) $(subst .cc,.hh,$<)+
	$(CXX)  -fPIC -O2  -shared $(CPPFLAGS) $< $(subst .cc,Dict.cc,$<) -o $@

clean:
	-rm $(TARGET) $(TARGET).o
	-rm $(addsuffix .o,$(MODULES))
	-rm *Dict.*
	-rm *.so
	-rm *\~

