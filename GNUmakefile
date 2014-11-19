# $Id: GNUmakefile,v 1.1.1.1 2008-06-04 13:59:29 mona Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := starliters
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
    G4INSTALL = /home/anna/geant4/geant4.9.6
endif


.PHONY: all
all: lib bin


#
# define G4LISTS_BASE, if you have your own physics lists area installed 
# point G4LISTS_BASE to the directory, that contains the subdirectory 'hadronic'.
#
ifndef G4LISTS_BASE
  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
  G4LISTS_BASE = $(G4INSTALL)/physics_lists
else
#  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
endif

#
# Select your physics lists to link against.
#
# EXTRALIBS += -lFTFC
# EXTRALIBS += -lFTFP
# EXTRALIBS += -lLBE
# EXTRALIBS += -lLHEP
# EXTRALIBS += -lLHEP_GN
# EXTRALIBS += -lLHEP_HP
# EXTRALIBS += -lLHEP_LEAD
# EXTRALIBS += -lLHEP_BERT_HP
# EXTRALIBS += -lLHEP_BIC_HP
# EXTRALIBS += -lLHEP_LEAD_HP
# EXTRALIBS += -lLHEP_PRECO
# EXTRALIBS += -lQGSP_BERT
#####EXTRALIBS += -lLHEP_PRECO_HP
# EXTRALIBS += -lQGSC
# EXTRALIBS += -lQGSC_LEAD
# EXTRALIBS += -lQGSC_LEAD_HP
# EXTRALIBS += -lQGSP
# EXTRALIBS += -lQGSP_GN
# EXTRALIBS += -lQGSP_HP
# EXTRALIBS += -lLHEP_BERT
# EXTRALIBS += -lLHEP_BIC
# EXTRALIBS += -lQGSP_BIC
#####EXTRALIBS += -lPackaging

##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/FTFC/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/FTFP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LBE/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BERT/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BIC/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_GN/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_HP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BERT_HP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BIC_HP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_LEAD/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_LEAD_HP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_PRECO/include
#####CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_PRECO_HP/include
#####CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/Packaging/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSC/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSC_LEAD/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSC_LEAD_HP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_BERT/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_BIC/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_GN/include
##CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_HP/include

##ROOT stuff:
CPPFLAGS += -I$(shell root-config --incdir)
EXTRALIBS = $(shell root-config --glibs)
CPPFLAGS += -g -ggdb -O0

#.PHONY: all
#all: lib bin

include $(G4INSTALL)/config/binmake.gmk
