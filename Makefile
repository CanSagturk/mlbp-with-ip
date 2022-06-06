
# If USE_CPP is set use this compiler, otherwise standard g++
ifdef USE_CPP
CPP=$(USE_CPP)
else
CPP=g++
endif

#GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
GIT_VERSION := "NO GIT"

# if CPLEXDIR is set, include CPLEX sources and link the libraries
ifdef CPLEXDIR
	CPLEX_FLAGS=-DCPLEX -I$(CPLEXDIR)/cplex/include -I$(CPLEXDIR)/concert/include -I$(CPLEXDIR)/cpoptimizer/include -DIL_STD
	CPLEX_LDFLAGS=-L$(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic -L$(CPLEXDIR)/concert/lib/x86-64_linux/static_pic -L$(CPLEXDIR)/cpoptimizer/lib/x86-64_linux/static_pic -lcp -lilocplex -lconcert -lcplex -lpthread -ldl
else
	CPLEX_FLAGS=-DNOCPLEX
	CPLEX_LDFLAGS=
endif

CCFLAGS=-pipe -Wall -DVERSION=\"$(GIT_VERSION)\" $(DEFS) -Wno-ignored-attributes -std=c++17

ifeq (,$(wildcard ./private.h))
CCFLAGS += -DNO_PRIVATE
endif

ifndef O4
CCFLAGS += -g
#CCFLAGS += -pg
else
CCFLAGS += -O3
CCFLAGS += -DNDEBUG
endif

.SUFFIXES: .o .cpp .d
.cpp.o: 
	$(CPP) -c $(CCFLAGS) $(CPLEX_FLAGS) $< -o $@
.cpp.d:
	$(CPP) -MM $(CCFLAGS) -MF $@ $<

HEADERS:=$(wildcard *.h) $(wildcard ./lib/*.h) $(wildcard ./mip/*.h)

SRCS=$(wildcard *.cpp) $(wildcard ./lib/*.cpp) $(wildcard ./mip/*.cpp)
DEPS=$(SRCS:.cpp=.d)
OBJS=$(SRCS:.cpp=.o)
DEPS=$(SRCS:.cpp=.d)

.PHONY: all clean

all: exes

mlbp: $(OBJS)
	$(CPP) $(CCFLAGS) -o mlbp $(OBJS) $(CPLEX_LDFLAGS) $(LDDEFS)

exes: mlbp 

clean:
	rm -f coating $(DEPS) $(OBJS)

sinclude $(DEPS)
