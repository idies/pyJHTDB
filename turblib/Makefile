# $Id: Makefile,v 1.9 2009-12-01 19:23:49 eric Exp $

#	Copyright 2011 Johns Hopkins University
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

OSARCH := $(shell uname -sp)

ifeq ($(OSARCH),Darwin i386)
# Compile both 32- and 64-bit code under MacOS X for Intel
# ARCH_FLAGS = -arch i386 -arch x86_64
else
	ARCH_FLAGS =
endif

RM     = rm -f
CFLAGS = -Wall
LDLIBS = -lm
CP     = cp
MKDIR  = mkdir -p

CC     = gcc -g $(ARCH_FLAGS)
FC     = gfortran $(ARCH_FLAGS)

OBJ =	soapC.o \
	soapClient.o \
	stdsoap2.o \
        turblib.o

all: DEMO_turbc DEMO_turbf DEMO_mhdc DEMO_mhdf DEMO_channelc DEMO_channelf DEMO_mixingc DEMO_mixingf DEMO_getCutoutc DEMO_getCutoutf

DEMO_getCutoutc : $(OBJ) DEMO_getCutoutc.o
	$(CC) -o $@ $(OBJ) DEMO_getCutoutc.o $(LDLIBS)

DEMO_getCutoutc.o : compiler_flags

DEMO_getCutoutf : $(OBJ) DEMO_getCutoutf.o
	 $(FC) -o $@ $(OBJ) DEMO_getCutoutf.o $(LDLIBS)

DEMO_getCutoutf.o : DEMO_getCutoutf.f90
	 $(FC) -c DEMO_getCutoutf.f90

DEMO_mhdc : $(OBJ) DEMO_mhdc.o
	 $(CC) -o $@ $(OBJ) DEMO_mhdc.o $(LDLIBS)

DEMO_mhdc.o: compiler_flags

DEMO_turbc : $(OBJ) DEMO_turbc.o
	 $(CC) -o $@ $(OBJ) DEMO_turbc.o $(LDLIBS)

DEMO_turbc.o: compiler_flags

DEMO_turbf : $(OBJ) DEMO_turbf.o
	 $(FC) -o $@ $(OBJ) DEMO_turbf.o $(LDLIBS)

DEMO_turbf.o : DEMO_turbf.f90
	 $(FC) -c DEMO_turbf.f90

DEMO_mhdf : $(OBJ) DEMO_mhdf.o
	 $(FC) -o $@ $(OBJ) DEMO_mhdf.o $(LDLIBS)

DEMO_mhdf.o : DEMO_mhdf.f90
	 $(FC) -c DEMO_mhdf.f90

DEMO_channelc : $(OBJ) DEMO_channelc.o
	 $(CC) -o $@ $(OBJ) DEMO_channelc.o $(LDLIBS)

DEMO_channelc.o: compiler_flags

DEMO_channelf : $(OBJ) DEMO_channelf.o
	 $(FC) -o $@ $(OBJ) DEMO_channelf.o $(LDLIBS)

DEMO_channelf.o : DEMO_channelf.f90
	 $(FC) -c DEMO_channelf.f90

DEMO_mixingc : $(OBJ) DEMO_mixingc.o
	 $(CC) -o $@ $(OBJ) DEMO_mixingc.o $(LDLIBS)

DEMO_mixingc.o: compiler_flags

DEMO_mixingf : $(OBJ) DEMO_mixingf.o
	 $(FC) -o $@ $(OBJ) DEMO_mixingf.o $(LDLIBS)

DEMO_mixingf.o : DEMO_mixingf.f90
	 $(FC) -c DEMO_mixingf.f90

stdsoap2.o: stdsoap2.c
	$(CC) $(CFLAGS) -c $<

static_lib: $(OBJ)
	ar rcs libJHTDB.a $(OBJ)

install: static_lib
	$(MKDIR) $(JHTDB_PREFIX)/include
	$(MKDIR) $(JHTDB_PREFIX)/lib
	$(CP) *.h $(JHTDB_PREFIX)/include/
	$(CP) libJHTDB.a $(JHTDB_PREFIX)/lib/

# Regenerate the gSOAP interfaces if required
TurbulenceService.h : wsdl

# Update the WSDL and gSOAP interfaces
wsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

testwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://test.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

mhdtestwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://mhdtest.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

devwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://dev.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

mhddevwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://mhddev.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

prodtestwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://prodtest.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

clean:
	$(RM) *.o *.exe turbf turbc mhdc mhdf channelc channelf mixingc mixingf compiler_flags getCutoutc getCutoutf

spotless: clean
	$(RM) soapClient.c TurbulenceServiceSoap.nsmap soapH.h TurbulenceServiceSoap12.nsmap soapStub.h soapC.c TurbulenceService.h

.SUFFIXES: .o .c .x

.c.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: force
compiler_flags: force
	echo '$(CFLAGS)' | cmp -s - $@ || echo '$(CFLAGS)' > $@

$(OBJ): compiler_flags

