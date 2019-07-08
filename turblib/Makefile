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

all: turbc turbf mhdc mhdf channelc channelf mixingc mixingf getCutoutc getCutoutf

getCutoutc : $(OBJ) getCutoutc.o
	$(CC) -o $@ $(OBJ) getCutoutc.o $(LDLIBS)

getCutoutc.o : compiler_flags

getCutoutf : $(OBJ) getCutoutf.o
	 $(FC) -o $@ $(OBJ) getCutoutf.o $(LDLIBS)

getCutoutf.o : getCutoutf.f90
	 $(FC) -c getCutoutf.f90

mhdc : $(OBJ) mhdc.o
	 $(CC) -o $@ $(OBJ) mhdc.o $(LDLIBS)

mhdc.o: compiler_flags

turbc : $(OBJ) turbc.o
	 $(CC) -o $@ $(OBJ) turbc.o $(LDLIBS)

turbc.o: compiler_flags

turbf : $(OBJ) turbf.o
	 $(FC) -o $@ $(OBJ) turbf.o $(LDLIBS)

turbf.o : turbf.f90
	 $(FC) -c turbf.f90

mhdf : $(OBJ) mhdf.o
	 $(FC) -o $@ $(OBJ) mhdf.o $(LDLIBS)

mhdf.o : mhdf.f90
	 $(FC) -c mhdf.f90

channelc : $(OBJ) channelc.o
	 $(CC) -o $@ $(OBJ) channelc.o $(LDLIBS)

channelc.o: compiler_flags

channelf : $(OBJ) channelf.o
	 $(FC) -o $@ $(OBJ) channelf.o $(LDLIBS)

channelf.o : channelf.f90
	 $(FC) -c channelf.f90

mixingc : $(OBJ) mixingc.o
	 $(CC) -o $@ $(OBJ) mixingc.o $(LDLIBS)

mixingc.o: compiler_flags

mixingf : $(OBJ) mixingf.o
	 $(FC) -o $@ $(OBJ) mixingf.o $(LDLIBS)

mixingf.o : mixingf.f90
	 $(FC) -c mixingf.f90

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

