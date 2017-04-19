#
# make : 18-04-2017 
#
# Description: Makefile for compiling IRK-Fixed-Point project.
# System: Ubuntu 16.04 lts (64 bits)
#
#


# Compiler.
SYS = Linux
CC = /usr/bin/gcc

# Wolfram Mathematica. 
MATHLINKPATH = $(MATHEMATICA_HOME)/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
MPREP = $(MATHLINKPATH)/mprep

# Compiler options.
LIBFLAGS=-L$(MATHLINKPATH)

FLAGS=-I$(MATHLINKPATH) -O2 -Wall -std=c99 -fno-common -mfma
OMP_FLAGS=-fopenmp

# Objects
OBJECTS = Terminal-IRKFP.o  Common-IRKFP.o Problems.o GaussCoefficients.o 
REALOBJECTS = math-IRKFPtm.o math-IRKFP.o Common-IRKFP.o Problems.o GaussCoefficients.o 

.PHONY: clean help 

all : 
	make clean;
	make term-IRKFP;
	make math-IRKFP;

term-IRKFP: $(OBJECTS) 
	${CC} $(OBJECTS) ${OMP_FLAGS} -lm -o Terminal-IRKFP.exe 
     
math-IRKFPtm.c : math-IRKFP.tm
	${MPREP} math-IRKFP.tm -o math-IRKFPtm.c

%.o: %.c def.h 
	${CC} ${FLAGS} ${OMP_FLAGS} -c -I./  $<				   

math-IRKFP : $(REALOBJECTS)    
	${CC}  $(REALOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -o math-IRKFP.exe 


clean :
	rm -rf *.o *~ core *tm.* *.exe 

help:
	@echo "Valid targets;"
	@echo " term-IRKFP"
	@echo " math-IRKFP"




