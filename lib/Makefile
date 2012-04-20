CXX = g++
CC  = gcc
FC  = gfortran
LD  = g++

# FC  = gfortran-mp-4.6

.SUFFIXES: .o .f90 .cpp .cu .s

OFLAGS = -O3 -g -m64 -Wall -msse3 -funroll-all-loops
CFLAGS = $(OFLAGS) -D__$(MACOSX)__ 
CXXFLAGS = $(CFLAGS) 
FCFLAGS = $(OFLAGS)

GXXLIB = -L/opt/local/lib/gcc46 -L/usr/lib  -lstdc++


OFLAGS += -D_FPESIG_ENABLE_ 

PROG = nbody_c
OBJS = nbody_c.o

PROG1 = nbody_f90
OBJS1 = nbody_f.o
OBJS1 = nbody_adaptive.o

LIBOBJ = libdh.o
LIBDH  = DH

# LIBS = -L/opt/local/lib -lCGAL  -lgmp
# LIBS += -lefence


all:  $(LIBDH) $(PROG) $(PROG1)

$(PROG): $(OBJS)
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) -L./ -l$(LIBDH)

$(PROG1): $(OBJS1)
	$(FC) $(OFLAGS) -o $@ $(OBJS1) -L./ -l$(LIBDH) $(GXXLIB)

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
