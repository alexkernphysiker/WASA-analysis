
#
# This can be used as a template Makefile for analysis programs
# which uses the RootSorter library
#
#
# This is where you have installed the RootSorter 

BASEDIR = $ROOTSORTERSYS/core

#
# TARGET is the basename of the file containing main()
# This should be a separate file 

TARGET   = main

#
# MODULES are all .cc/.hh file combinations containing your
# own classes
# The ROOT linkage has to be specified in 'LinkDef.hh'
#

MODULES  =  examplejob

# 
# Starting from here no changes should be necessary 
# -lWasaRecFD 

CPPFLAGS := `$(BASEDIR)/bin/sorter-config -cpp` -g  
CXXFLAGS := `$(BASEDIR)/bin/sorter-config -cxx` -fno-strict-aliasing  -g 
LDFLAGS  := `$(BASEDIR)/bin/sorter-config -ld -libs-wasa -libs-wasa-ana` -lWasaRecFD -lWasaRecSE -lWasaRecPS -lWasaRecMDC -lWasaRecCD -lWasaParameter -lWasaAnaRaw -lWasaRecFPC -lvt -g 

all: $(TARGET)

$(TARGET): $(TARGET).o 
	$(CXX) $(LDFLAGS) $^ -o $@

ifneq (xx,x$(MODULES)x)
$(TARGET): $(TARGET).so

$(TARGET).so: $(addsuffix .o,$(MODULES)) $(TARGET)Dict.o
	$(CXX) -fPIC -O2 -shared $^ -o $@

$(TARGET)Dict.cc: $(addsuffix .hh,$(MODULES)) LinkDef.hh
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



