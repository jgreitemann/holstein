DEFINES+= -DMCL_DUMP_BUFFER=0
DEFINES+= -DMCL_MEASUREMENTS_APPEND=0
DEFINES+= -DMCL_MCL_RNG_MT
DEFINES+= -DMEASURE_KIN_ENERGY

OBJS_MPI = dump.mpi.o parser.mpi.o measurements.mpi.o evalable.mpi.o observable.mpi.o random.mpi.o mc.mpi.o main.mpi.o runner.mpi.o
OBJS_SINGLE = dump.single.o parser.single.o measurements.single.o evalable.single.o observable.single.o random.single.o mc.single.o main.single.o runner_single.single.o
DEFINES_SINGLE = $(DEFINES) -DMCL_SINGLE
OBJS_PT = dump.pt.o parser.pt.o measurements.pt.o evalable.pt.o observable.pt.o random.pt.o mc.pt.o main.pt.o runner.pt.o
DEFINES_PT = -DMCL_PT

MCLL  = $(MC_CODE_DIR)/load_leveller/trunk
APPMCLL = $(MC_CODE_DIR)/holstein

CC_MPI = $(MPICC)
LD_MPI = $(MPICC)
CC_SINGLE = g++
LD_SINGLE = g++
CC_PT = $(MPICC)
LD_PT = $(MPICC)
CFLAGS = -O3 -Wno-deprecated --short-enums -g -ansi -ffast-math -Wall
CFLAGS_MPI  = $(CFLAGS) $(DEFINES)
CFLAGS_SINGLE  = $(CFLAGS) $(DEFINES_SINGLE)
CFLAGS_PT  = $(CFLAGS) $(DEFINES_PT)

INCLUDE = -I$(MCLL) -I$(APPMCLL) 
LDFLAGS = -lm

RM = /bin/rm -f

all: mc mc_single cleano

mc : $(OBJS_MPI)
	$(LD_MPI) $(LDFLAGS) -o $@ $(OBJS_MPI)

%.mpi.o : $(APPMCLL)/%.cpp
	$(CC_MPI) -c $(CFLAGS_MPI) $(INCLUDE) $< -o $@

%.mpi.o : $(MCLL)/%.cpp
	$(CC_MPI) -c $(CFLAGS_MPI) $(INCLUDE) $< -o $@

mc_single : $(OBJS_SINGLE)
	$(LD_SINGLE) $(LDFLAGS) -o $@ $(OBJS_SINGLE)

%.single.o : $(APPMCLL)/%.cpp
	$(CC_SINGLE) -c $(CFLAGS_SINGLE) $(INCLUDE) $< -o $@

%.single.o : $(MCLL)/%.cpp
	$(CC_SINGLE) -c $(CFLAGS_SINGLE) $(INCLUDE) $< -o $@

mc_pt : $(OBJS_PT)
	$(LD_PT) $(LDFLAGS) -o $@ $(OBJS_PT)

%.pt.o : $(APPMCLL)/%.cpp
	$(CC_PT) -c $(CFLAGS_PT) $(INCLUDE) $< -o $@

%.pt.o : $(MCLL)/%.cpp
	$(CC_PT) -c $(CFLAGS_PT) $(INCLUDE) $< -o $@

clean:
	$(RM) mc mc_single mc_pt $(OBJS_MPI) $(OBJS_SINGLE) $(OBJS_PT)

cleano:
	$(RM) $(OBJS_MPI) $(OBJS_SINGLE) $(OBJS_PT)



