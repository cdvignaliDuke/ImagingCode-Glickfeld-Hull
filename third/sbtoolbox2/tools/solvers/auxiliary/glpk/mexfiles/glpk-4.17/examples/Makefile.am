# Process this file with automake to produce Makefile.in

INCLUDES = -I../include

LDADD = ../src/libglpk.la

bin_PROGRAMS = glpsol tspsol

glpsol_SOURCES = glpsol.c

tspsol_SOURCES = tspsol.c

check: glpsol$(EXEEXT)
	./glpsol$(EXEEXT) --mps ./plan.mps
