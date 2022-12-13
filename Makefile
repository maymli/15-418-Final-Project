OUTPUTDIR := bin/

CFLAGS := -std=c++14 -fvisibility=hidden -lpthread -Wall

ifeq (,$(CONFIGURATION))
	CONFIGURATION := release
endif

ifeq (debug,$(CONFIGURATION))
CFLAGS += -g
else ifeq (MPI, $(CONFIGURATION))
CFLAGS += -O2
else
CFLAGS += -O2 -fopenmp
endif

ifeq (MPI, $(CONFIGURATION))
SOURCES = openmpi-coloringv2.cpp
else
SOURCES := src/*.cpp
endif
HEADERS := src/*.h

ifeq (MPI, $(CONFIGURATION))
CXX = mpic++
endif
 
TARGETBIN := color-$(CONFIGURATION)

.SUFFIXES:
.PHONY: all clean

all: $(TARGETBIN)

$(TARGETBIN): $(SOURCES) $(HEADERS)
	$(CXX) -o $@ $(CFLAGS) $(SOURCES)

format:
	clang-format -i src/*.cpp src/*.h

clean:
	rm -rf ./color-$(CONFIGURATION)

# add a check in later
