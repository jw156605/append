# ***********************************************************************
#
#  Makefile for AppEnD
#  =====================================
#  Josh Welch
#  6/9/2014
# ***********************************************************************
AppEnD:AppEnD.o
	g++ -L./bamtools/lib AppEnD.o -o AppEnD -Wl,-rpath,./bamtools/lib -lbamtools
AppEnD.o:AppEnD.cc AppEnD.h
	mkdir ./bamtools/build; cd ./bamtools/build; cmake ..; make
	g++ -c -O3 AppEnD.cc -I./bamtools/include