# $Id: GNUmakefile,v 1.2 2014-11-25 11:52:10 Abhijit  Exp $
# --------------------------------------------------------------
# --------------------------------------------------------------

name := vegaPB01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

