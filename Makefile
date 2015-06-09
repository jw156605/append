# ***********************************************************************
#
#  Makefile for AppEnD
#  =====================================
#  Josh Welch
#  11/25/2014
# ***********************************************************************
AppEnD:AppEnD.o
	g++ -L/proj/prins_lab/Josh/pezmaster31-bamtools-9527e22/lib AppEnD.o -o AppEnD -Wl,-rpath,/proj/prins_lab/Josh/pezmaster31-bamtools-9527e22/lib -lbamtools
AppEnD.o:AppEnD.cc AppEnD.h
	g++ -c -O3 AppEnD.cc -I/proj/prins_lab/Josh/pezmaster31-bamtools-9527e22/include
