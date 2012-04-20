CXX = g++
CC  = gcc
FC  = gfortran
LD  = g++

# CXX = g++-mp-4.6
# CC  = gcc-mp-4.6
# FC  = gfortran-mp-4.6
# LD  = g++-mp-4.6

.SUFFIXES: .o .f90 .cpp .cu .s

OFLAGS = -O3 -g -m64 -Wall -msse3 -funroll-all-loops
# OFLAGS = -O0 -g -m64 -Wall -msse3 
# OFLAGS += -fopenmp
CFLAGS = $(OFLAGS) -I$(CUDA_TK)/include  -D__$(MACOSX)__  -I/opt/local/include
# CFLAGS += -D_CH2SIG_
CXXFLAGS = $(CFLAGS) 

GFORTRANLIBPATH = /opt/local/lib/gcc46


OFLAGS += -D_FPESIG_ENABLE_ 

PROG = nbody_c
OBJS = nbody_c.o

PROG1 = nbody_f90
OBJS1 = nbody_f.o

LIBOBJ = libdh.o
LIBDH  = DH

# LIBS = -L/opt/local/lib -lCGAL  -lgmp
# LIBS += -lefence


all:  $(LIBDH) $(PROG) $(PROG1)

$(PROG): $(OBJS)
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) -L./ -l$(LIBDH)

$(PROG1): $(OBJS1)
	$(LD) $(CXXFLAGS) -o $@ $(OBJS1) -L./ -l$(LIBDH) -L$(GFORTRANLIBPATH) -lgfortran

$(LIBDH): $(LIBOBJ)
	/bin/rm -rf lib$@.a
	ar qv lib$@.a $^
	ranlib lib$@.a


.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

.f90.o: 
	$(FC) $(FCFLAGS) -c $< -o $@

.cpp.s: 
	$(CXX) $(CXXFLAGS) -S $< -o $@

clean:
	/bin/rm -rf *.o version.h
	/bin/rm -rf $(PROG) $(PROG1) lib$(VLIB).a

$(OBJS) : libdh.h libdh_inc.h particle.h vector3.h kepler.h dh.h  mytimer.h
