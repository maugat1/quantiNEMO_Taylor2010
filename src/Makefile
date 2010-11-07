################################################################################
# Makefile for quantiNEMO
################################################################################

# whether this is a debugging build or not
#DEBUG		:= YES
DEBUG		:= NO

################################################################################
# Compiler and linker
################################################################################

CC   	:= g++
CXX    	:= g++

################################################################################
# Preprocessor directives
################################################################################

DEBUG_DEFS	:= -D_DEBUG
RELEASE_DEFS	:=

################################################################################
# Include paths, library flags, and libraries
################################################################################

INCL		:= -I.
LDFLAGS		:=
LDLIBS		:=

################################################################################
# Compiler flags
################################################################################

DEBUG_CXXFLAGS     := -Wall -g -fPIC
RELEASE_CXXFLAGS   := -Wall -O3 -fPIC

ifeq (YES, ${DEBUG})
    BASE_CXXFLAGS       := ${DEBUG_CXXFLAGS} ${DEBUG_DEFS}
else
    BASE_CXXFLAGS       := ${RELEASE_CXXFLAGS} ${RELEASE_DEFS}
endif

CXXFLAGS := ${BASE_CXXFLAGS} ${INCL}

################################################################################
# Source files
################################################################################

SRCS := $(wildcard *.cpp) $(wildcard *.cc)

################################################################################
# Targets of the build
################################################################################

OBJS 	:= $(sort $(addsuffix .o, $(basename ${SRCS})))
BIN 	:= quantiNemo
LIB	:= $(addprefix lib, $(addsuffix .so, ${BIN}))

################################################################################
# Rules
################################################################################

all: bin

# compiles the objects (default rules will take care of the .o files)
objects: ${OBJS}

# compiles the executable
bin: ${BIN}

${BIN}: objects
	${CXX} ${OBJS} -o $@ ${CXXFLAGS} ${LDFLAGS} ${LDLIBS}

# compiles a shared library 
lib: ${LIB}

${LIB}: objects
	${CXX} -shared ${OBJS} -o $@ ${CXXFLAGS} ${LDFLAGS} ${LDLIBS}

# cleans up created files
clean:
	-rm -f ${OBJS} ${BIN} ${LIB}

# determines dependencies
depend:
	${CC} -M ${SRCS} > $@ ${COMDIR}
#	makedepend ${INCL} ${SRCS}	

#include depend
