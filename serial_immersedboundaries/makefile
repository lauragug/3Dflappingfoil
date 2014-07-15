#SWI= -c -r8 -arch ev6 -fast -tune ev6 -O3 (f90xime054)
#SWI= -c -O3 -w -Vaxlib #-tpp7 -xM
#SWI= -c -mcmodel=large -O3 -w -Vaxlib #-tpp7 -xM
SWI= -c  -O3 -w -Vaxlib #-tpp7 -xMcd
cogetn2: cogenn.o cogetn.o cogehd.o cogebc-rv.o forces.o cogefopi.o pressnu.o ffts.o fishpack.o cogeso.o lesmodc.o  mmovie.o stamod.o
	ifort -Vaxlib -o cogetn2  cogenn.o cogetn.o cogehd.o cogebc-rv.o forces.o cogefopi.o pressnu.o ffts.o fishpack.o cogeso.o lesmodc.o mmovie.o stamod.o
cogetn.o: cogetn.f param.f
	ifort $(SWI) cogetn.f
cogehd.o: cogehd.f param.f
	ifort $(SWI) cogehd.f
cogenn.o: cogenn.f param.f
	ifort $(SWI) cogenn.f
cogebc-rv.o: cogebc-rv.f param.f
	ifort $(SWI) cogebc-rv.f
forces.o: forces.f param.f
	ifort $(SWI) forces.f
cogefopi.o: cogefopi.f param.f
	ifort $(SWI) cogefopi.f
pressnu.o: pressnu.f param.f
	ifort $(SWI) pressnu.f
ffts.o: ffts.f
	ifort $(SWI) ffts.f
fishpack.o: fishpack.f
	ifort $(SWI) fishpack.f
cogeso.o: cogeso.f param.f
	ifort $(SWI) cogeso.f
lesmodc.o: lesmodc.f param.f
	ifort $(SWI) lesmodc.f
mmovie.o: mmovie.f param.f
	ifort $(SWI) mmovie.f
stamod.o: stamod.f param.f
	ifort $(SWI) stamod.f
clean:
	rm -f *.o
