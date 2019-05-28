# ***********************************************************************
#
#  Makefile for AppEnD
#  =====================================
#  Josh Welch
#  6/9/2014
# ***********************************************************************
AppEnD:AppEnD.o
	g++ -L./bamtools/build/src/api AppEnD.o -o AppEnD -Wl,-rpath,./bamtools/lib -lbamtools -lz
AppEnD.o:AppEnD.cc AppEnD.h
	mkdir ./bamtools/build; cd ./bamtools/build; cmake ..; make
	g++ -c -O3 AppEnD.cc -I./bamtools/src
