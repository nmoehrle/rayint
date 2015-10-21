MVE_ROOT ?= ../mve
TARGET := all
include ${MVE_ROOT}/Makefile.inc

CXXFLAGS += -I${MVE_ROOT}/libs -std=c++11 -fopenmp
LDLIBS += -lpng -ltiff -ljpeg -fopenmp

COMMON := libmve.a libmve_util.a

all: raycast
raycast: raycast.o $(COMMON)

clean:
	${RM} raycast *.o Makefile.dep

.PHONY: clean
