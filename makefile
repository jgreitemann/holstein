DEFINES+= -DMCL_DUMP_BUFFER=0
DEFINES+= -DMCL_MEASUREMENTS_APPEND=0
DEFINES+= -DMCL_MCL_RNG_MT

MODE=MPI
#MODE=SINGLE
#MODE=PT

OBJS = dump.o parser.o measurements.o evalable.o observable.o random.o mc.o main.o
OBJSLN = dump.LN.o parser.LN.o measurements.LN.o evalable.LN.o observable.LN.o random.LN.o mc.LN.o runner_single.LN.o merge.LN.o

ifeq ($(MODE),MPI)
  OBJS+=runner.o
endif

ifeq ($(MODE),SINGLE)
  OBJS+=runner_single.o
  DEFINES+= -DMCL_SINGLE
endif

ifeq ($(MODE),PT)
  OBJS+=runner.o
  DEFINES+= -DMCL_PT
endif

MCLL  = $(MC_CODE_DIR)/load_leveller/trunk
APPMCLL = $(MC_CODE_DIR)/holstein

CC = $(MPI)/bin/mpiCC
LD = $(MPI)/bin/mpiCC
ifeq ($(MODE),SINGLE)
  CC=g++
  LD=g++
endif
CFLAGS  = -O3 -Wno-deprecated -g -ansi -ffast-math -Wall $(DEFINES)
INCLUDE = -I$(MCLL) -I$(APPMCLL) 
LDFLAGS = -lm  
SUPERLP = 

CCLN = g++
LDLN = g++
CFLAGSLN  = $(CFLAGS) -DMCL_SINGLE
INCLUDELN = $(INCLUDE)
LDFLAGSLN = $(LDFLAGS)
SUPERLPLN = $(SUPERLP)

RM = /bin/rm -f

#all: mc merge cleano
all: mc cleano

mc : $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(SUPERLP)

merge : $(OBJSLN)
	$(LDLN) $(LDFLAGSLN) -o $@ $(OBJSLN) $(SUPERLPLN)

%.o : $(APPMCLL)/%.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

%.LN.o : $(APPMCLL)/%.cpp
	$(CCLN) -c $(CFLAGSLN) $(INCLUDELN)  $< -o $@

%.o : $(MCLL)/%.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

%.LN.o : $(MCLL)/%.cpp
	$(CCLN) -c $(CFLAGSLN) $(INCLUDELN) $< -o $@ 

clean:
	$(RM) mc merge $(OBJS) $(OBJSLN)

cleano:
	$(RM) $(OBJS) $(OBJSLN)



