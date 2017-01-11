#
# make : 01-10-2016 
#
# Description: Makefile for compiling IRK-Fixed-Point project.
#


# Compiler.
SYS = Linux
CC = /usr/bin/gcc

# Wolfram Mathematica. 
MATHLINKPATH = $(MATHEMATICA_HOME)/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
MPREP = $(MATHLINKPATH)/mprep

# Compiler options.
LIBFLAGS=-L$(MATHLINKPATH)

ifeq ($(FMA),NO)
  FLAGS=-I$(MATHLINKPATH) -O2 -Wall -std=c99 -fno-common
else
  FLAGS=-I$(MATHLINKPATH) -O2 -Wall -std=c99 -fno-common -mfma 
endif  

FLAGS=-I$(MATHLINKPATH) -O2 -Wall -std=c99 -fno-common
OMP_FLAGS=-fopenmp

# Objects
OBJECTS = GaussTerminal.o GaussInitData.o GaussCommon.o GaussUserProblem.o GaussCoefficients.o 
REALOBJECTS = math-Gausstm.o math-Gauss.o GaussCommon.o GaussUserProblem.o GaussCoefficients.o 
QUADOBJECTS = quad-Gausstm.o quad-Gauss.o GaussCommon.o GaussUserProblem.o GaussCoefficients.o 
OBJECTS3 = math-GaussParatm.o math-GaussPar.o GaussCommon.o GaussUserProblem.o  GaussCoefficients.o 


.PHONY: clean help 

all : 
	make clean;
	make term-Gauss;
	make math-Gauss;
	make quad-Gauss;
	make math-GaussPar

term-Gauss: $(OBJECTS) 
	${CC} $(OBJECTS) ${OMP_FLAGS} -lm -lquadmath -o GaussTerminal.exe      
math-Gausstm.c : math-Gauss.tm
	${MPREP} math-Gauss.tm -o math-Gausstm.c

quad-Gausstm.c : quad-Gauss.tm
	${MPREP} quad-Gauss.tm -o quad-Gausstm.c

math-GaussParatm.c : math-GaussPar.tm
	${MPREP} math-GaussPar.tm -o math-GaussParatm.c

%.o: %.c prec.h def.h 
	${CC} ${FLAGS} ${OMP_FLAGS} -c -I./  $<				   

math-Gauss : $(REALOBJECTS)    
	${CC}  $(REALOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -o math-Gauss.exe 

quad-Gauss : $(QUADOBJECTS)   
	${CC} $(QUADOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -lquadmath -o quad-Gauss.exe

math-GaussPar : $(OBJECTS3) 
	${CC} $(OBJECTS3) $(LIBFLAGS) ${OMP_FLAGS}  -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -lquadmath -o math-GaussPar.exe

clean :
	rm -rf *.o *~ core *tm.* *.exe 

help:
	@echo "Valid targets;"
	@echo " term-Gauss"
	@echo " math-Gauss"
	@echo " quad-Gauss"
	@echo " math-GaussPar"



