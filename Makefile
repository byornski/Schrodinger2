FC:=gfortran
DEBUG:=0
FFLAGS:=-O1 -g -fcheck=all  -fbacktrace -std=f2008 -DDEBUG=${DEBUG}
WARNINGS=-pedantic -Wall -Wextra -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface -Wunused-parameter -ffpe-trap=invalid,zero,overflow
BLD_DIR=build
SRC_DIR=src
MOD_DIR=build

COMPILE=${FC} ${FFLAGS} ${WARNINGS} -I${MOD_DIR} -J${MOD_DIR}


all: ham.exe

%.exe: ${BLD_DIR}/%.o
	${COMPILE} -o $@ ${BLD_DIR}/*.o

${BLD_DIR}/%.o: ${SRC_DIR}/%.f90 
	${COMPILE} -c -o $@ $<

#${BLD_DIR}/shooting.o: ${BLD_DIR}/rungekutta.o ${BLD_DIR}/well_parameters.o
#${BLD_DIR}/rungekutta.o: ${BLD_DIR}/well_parameters.o

${BLD_DIR}/ham.o: ${BLD_DIR}/wave.o ${BLD_DIR}/util.o ${BLD_DIR}/hamiltonian.o
${BLD_DIR}/hamiltonian.o: ${BLD_DIR}/wave.o
${BLD_DIR}/potential.o: ${BLD_DIR}/well_parameters.o ${BLD_DIR}/basis.o ${BLD_DIR}/util.o
${BLD_DIR}/wave.o: ${BLD_DIR}/basis.o ${BLD_DIR}/potential.o ${BLD_DIR}/util.o
${BLD_DIR}/basis.o: ${BLD_DIR}/util.o


.PHONY: clean
clean:
	rm -f build/* ham.exe
