#SWI= -c  -O3 -w -Vaxlib #-tpp7 -xM
SWI= -g -r8
SHELL = /bin/sh
F90   = /usr/rels/mpich-1.2.5..9/bin/mpif90
MAKE  = make --no-print-directory
MPIRUN= mpirun
NP_ARG= -np 
#includedir=/usr/rels/mpich-1.2.5..9/include

all: all-redirect

PROGRAMS = prova

prova: cogenn.f outfile.f cogebc-rv.f cogefopi.f cogehd.f cogeso.f cogetn.f exchng.f ffts.f fishpack.f pressnu.f90 solq.f90 forces.f
	${F90} ${SWI} -o prova cogenn.f outfile.f cogebc-rv.f cogefopi.f cogehd.f cogeso.f cogetn.f exchng.f ffts.f fishpack.f pressnu.f90 solq.f90 forces.f

all-redirect: ${PROGRAMS}

clean:
	rm -f *.o

run: ${PROGRAMS}
	${MPIRUN} ${NP_ARG} 4 prova
