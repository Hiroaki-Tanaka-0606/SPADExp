# CC      = pgcc -Mbounds -Ktrap=fp
# CC      = gcc -std=c99 -O3 
# CFLAGS  = -g
#CC      = icc -openmp -O2 -I/usr/local/include -I/home/ozaki/include
#CC      = gcc -std=c99 -O3 -I/usr/local/include -I/home/ozaki/include
#CC      = gcc -std=c99 -O3 -I/usr/local/include 
#LIB     = -llapack -lblas

MKLROOT = /opt/intel/mkl/lib/intel64
CC      = gcc -Dnoomp -std=c99 -O3 -I/usr/local/include -mcmodel=large
LIB     =  -Wl,--start-group ${MKLROOT}/libmkl_intel_lp64.a ${MKLROOT}/libmkl_intel_thread.a ${MKLROOT}/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -ldl -lm 

CFLAGS  = -g
LDFLAGS = 

OBJS    = adpack.o Inputtools.o readfile.o Set_Init.o addfunc.o \
          All_Electron.o Total_Energy.o \
          Density.o Density_V.o Density_PCC.o Core.o Hartree.o \
          HokanF.o XC_Xa.o XC_CA.o XC_VWN.o \
          XC4atom_PBE.o VP.o Hamming_O.o Hamming_I.o XC_PBE.o \
          XC_PW91C.o XC_EX.o Gauss_LEQ.o Gauss_Legendre.o BHS.o \
          TM.o MBK.o MR.o PAO_RadialF.o MPAO_RadialF.o \
          Frho_V.o VNLF.o Output.o Initial_Density.o GR_Pulay.o \
          Simple_Mixing.o Multiple_PAO.o Calc_Vlocal.o Init_VPS.o \
          Generate_VNL.o Log_DeriF.o ghost.o Restart.o Make_EDPP.o \
          Find_LESP.o Empty_VPS.o QuickSort.o mimic_omp.o \
          FEM_All_Electron.o FEMLDA_All_Electron.o \
          FEMHF_All_Electron.o FEMHF_ERI.o FEMHF_JKLM.o \
          All_Electron_NSCF.o

PROG    = adpack
DESTDIR = ../work

all:	$(PROG)

adpack: $(OBJS)
	$(CC) $(OBJS) $(LIB) -o adpack $(LDFLAGS)
#
#
#
adpack.o: adpack.c adpack.h
	$(CC) -c adpack.c
Inputtools.o: Inputtools.c
	$(CC) -c Inputtools.c
readfile.o: readfile.c adpack.h Inputtools.h
	$(CC) -c readfile.c
Set_Init.o: Set_Init.c adpack.h
	$(CC) -c Set_Init.c
addfunc.o: addfunc.c adpack.h
	$(CC) -c addfunc.c
All_Electron.o: All_Electron.c adpack.h
	$(CC) -c All_Electron.c
Total_Energy.o: Total_Energy.c adpack.h
	$(CC) -c Total_Energy.c
Density.o: Density.c adpack.h
	$(CC) -c Density.c
Density_V.o: Density_V.c adpack.h
	$(CC) -c Density_V.c
Density_PCC.o: Density_PCC.c adpack.h
	$(CC) -c Density_PCC.c
Core.o: Core.c adpack.h
	$(CC) -c Core.c
HokanF.o: HokanF.c adpack.h
	$(CC) -c HokanF.c
Hartree.o: Hartree.c adpack.h
	$(CC) -c Hartree.c
XC_Xa.o: XC_Xa.c adpack.h
	$(CC) -c XC_Xa.c
XC_CA.o: XC_CA.c adpack.h
	$(CC) -c XC_CA.c
XC_VWN.o: XC_VWN.c adpack.h
	$(CC) -c XC_VWN.c
XC4atom_PBE.o: XC4atom_PBE.c adpack.h
	$(CC) -c XC4atom_PBE.c
VP.o: VP.c adpack.h
	$(CC) -c VP.c
Hamming_O.o: Hamming_O.c adpack.h
	$(CC) -c Hamming_O.c
Hamming_I.o: Hamming_I.c adpack.h
	$(CC) -c Hamming_I.c
XC_PBE.o: XC_PBE.c adpack.h
	$(CC) -c XC_PBE.c
XC_PW91C.o: XC_PW91C.c adpack.h
	$(CC) -c XC_PW91C.c
XC_EX.o: XC_EX.c adpack.h
	$(CC) -c XC_EX.c
Gauss_LEQ.o: Gauss_LEQ.c adpack.h
	$(CC) -c Gauss_LEQ.c
Gauss_Legendre.o: Gauss_Legendre.c adpack.h
	$(CC) -c Gauss_Legendre.c
BHS.o: BHS.c adpack.h
	$(CC) -c BHS.c
TM.o: TM.c adpack.h
	$(CC) -c TM.c
MBK.o: MBK.c adpack.h
	$(CC) -c MBK.c
MR.o: MR.c adpack.h
	$(CC) -c MR.c
PAO_RadialF.o: PAO_RadialF.c adpack.h
	$(CC) -c PAO_RadialF.c
MPAO_RadialF.o: MPAO_RadialF.c adpack.h
	$(CC) -c MPAO_RadialF.c
Frho_V.o: VNLF.c adpack.h
	$(CC) -c Frho_V.c
VNLF.o: VNLF.c adpack.h
	$(CC) -c VNLF.c
#E_NL.o: E_NL.c adpack.h
#	$(CC) -c E_NL.c
Output.o: Output.c adpack.h
	$(CC) -c Output.c
Initial_Density.o: Initial_Density.c adpack.h
	$(CC) -c Initial_Density.c
Simple_Mixing.o: Simple_Mixing.c adpack.h
	$(CC) -c Simple_Mixing.c
GR_Pulay.o: GR_Pulay.c adpack.h
	$(CC) -c GR_Pulay.c
Multiple_PAO.o: Multiple_PAO.c adpack.h
	$(CC) -c Multiple_PAO.c
Calc_Vlocal.o: Calc_Vlocal.c adpack.h
	$(CC) -c Calc_Vlocal.c
Init_VPS.o: Init_VPS.c adpack.h
	$(CC) -c Init_VPS.c
Generate_VNL.o: Generate_VNL.c adpack.h
	$(CC) -c Generate_VNL.c
Log_DeriF.o: Log_DeriF.c adpack.h
	$(CC) -c Log_DeriF.c
ghost.o: ghost.c adpack.h
	$(CC) -c ghost.c
Restart.o: Restart.c adpack.h
	$(CC) -c Restart.c
Make_EDPP.o: Make_EDPP.c adpack.h
	$(CC) -c Make_EDPP.c
Find_LESP.o: Find_LESP.c adpack.h
	$(CC) -c Find_LESP.c
Empty_VPS.o: Empty_VPS.c adpack.h
	$(CC) -c Empty_VPS.c
QuickSort.o: QuickSort.c adpack.h
	$(CC) -c QuickSort.c
mimic_omp.o: mimic_omp.c mimic_omp.h
	 $(CC) -c mimic_omp.c
#
#
FEM_All_Electron.o: FEM_All_Electron.c adpack.h
	$(CC) -c FEM_All_Electron.c
FEMLDA_All_Electron.o: FEMLDA_All_Electron.c adpack.h
	$(CC) -c FEMLDA_All_Electron.c
FEMHF_All_Electron.o: FEMHF_All_Electron.c adpack.h FEMHF_ERI.h
	$(CC) -c FEMHF_All_Electron.c 
FEMHF_ERI.o: FEMHF_ERI.c FEMHF_ERI.h FEMHF_JKLM.h
	$(CC) -c FEMHF_ERI.c
FEMHF_JKLM.o: FEMHF_JKLM.c FEMHF_JKLM.h
	$(CC) -c FEMHF_JKLM.c 

#
install:	$(PROG)
	strip $(PROG)
	cp $(PROG) $(DESTDIR)/$(PROG)

clean:	
	rm -f $(PROG) $(OBJS)
