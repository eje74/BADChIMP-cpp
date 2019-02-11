
VERSION = -std=c++11

PROGRAM = IO-test
OS := $(shell uname)
HOST := $(shell hostname)

#export OMPI_CXX=g++-8
MPICC = mpic++ # build parallel c++ code

### turn on MPI by default
MPI = 1


### GNU gcc compiler flags
CFLAGS = -O0 -march=native -Wall -g
LDFLAGS += -lm

ifeq ($(OS), Darwin)
LIBS = -I/usr/include/malloc -Iinclude
endif

ifeq ($(MPI),1)
CC = $(MPICC)
CFLAGS += -D_MPI_
endif

CFLAGS += $(VERSION)
LDFLAGS += $(VERSION)

SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.P)

$(PROGRAM):	$(OBJS)
		$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

%.o : %.cpp
	$(CC) -MMD -c -o $@ $< $(LIBS) $(CFLAGS) 
	@cp $*.d $*.P; \
	  sed -e 's/\#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	  rm -f $*.d

### print variable values by command: make print-VARNAME
print-%  : ; @echo $* = $($*)

### copy binary to case-specific folder
flow grad_P perc source calc_magn_replace two-phase:
	mkdir -p $(DESTDIR)/$@/ && cp $(PROGRAM) $(DESTDIR)/$@/

.PHONY: clean
clean:
	@-rm -f src/*.o src/*.P src/*.d

### avoid including DEPS during clean (why create just to delete them?)
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS) 
endif

# end of Makefile 

