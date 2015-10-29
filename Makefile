.EXPORT_ALL_VARIABLES:

.PHONY: clean all

LIB_DIR = $(HOME)/lib

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

COMMON_DIR = $(HOME)/common


BASELIBS  = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR)
ALLIBS  = $(BASELIBS) -lCommandLineInterface -lKinematics

CPP             = g++
#CFLAGS	        = -g -O3 $(ROOTCFLAGS)
CFLAGS		= -pedantic -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC
DFLAGS		= -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC

INCLUDES        = -I./ -I$(COMMON_DIR) 
LFLAGS		= -g -fPIC
LIBS 		= $(ALLIBS)

CFLAGS += -Wl,--no-as-needed
LFLAGS += -Wl,--no-as-needed
DFLAGS += -Wl,--no-as-needed

O_FILES = Nucleus.o \
	Kinematics.o 

LIBRARIES = $(LIB_DIR)/libCommandLineInterface.so 

all: Reaction
	echo Done

Reaction: Reaction.cc $(O_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) $^ $(LIBS) -o $@
	cp Reaction $(HOME)/bin

lib%.so: Nucleus.o NucleusDictionary.o Kinematics.o KinematicsDictionary.o
	$(CPP) $(LFLAGS) -shared -Wl,-soname,$@ -o $(LIB_DIR)/$@ $^ $(BASELIBS) -lc

%.o: %.cc %.hh
	@echo Default .o rule
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

NucleusDictionary.o: NucleusDictionary.cc NucleusDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

NucleusDictionary.cc: Nucleus.hh NucleusLinkDef.h 
	 rm -f NucleusDictionary.cc NucleusDictionary.h; rootcint -f $@ -c Nucleus.hh NucleusLinkDef.h 

KinematicsDictionary.o: KinematicsDictionary.cc KinematicsDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

KinematicsDictionary.cc: Kinematics.hh KinematicsLinkDef.h 
	 rm -f KinematicsDictionary.cc KinematicsDictionary.h; rootcint -f $@ -c Kinematics.hh KinematicsLinkDef.h 

clean:
	rm *.o Reaction

tar:
	@echo "creating zipped tar-ball ... "
	tar -chvzf Reaction.tar.gz Makefile mass.dat \
	*.hh *.cc
