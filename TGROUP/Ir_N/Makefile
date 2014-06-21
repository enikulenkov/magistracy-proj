BUILD_HOME=../..
include $(BUILD_HOME)/Options

#PAGMO_PATH=/home/deimos/src/tests/genetic_algorithm/PaGMO/pagmo-code
FORTRAN_CCFLAGS = -g -fdollar-ok

BIN1 = sc_qdm_n
BIN2 = ga
BINS = $(BIN1) $(BIN2)

all: $(BIN1) $(BIN2)

SRCS0 = ga_utils_f.cpp
OBJS0 = ${SRCS0:.cpp=.o}

$(BIN1): SC_QDM_N.f90 $(OBJS0) 
	$(FORTRAN_CC) -o $@ $(FORTRAN_CCFLAGS) $+ ../../lib/libutils.a -lstdc++ -lopenbabel 

#INCPATH+=-I/usr/include/openbabel-2.0 -I$(PAGMO_PATH)/src
INCPATH+=-I/usr/include/openbabel-2.0

SRCS2 = ga.cpp sutton_chen_pot.cpp ga_utils.cpp
OBJS2 = ${SRCS2:.cpp=.o}

#LDFLAGS_USER =-L$(PAGMO_PATH)/src -lpagmo -lboost_serialization -lboost_system -lboost_thread -lboost_python -lpthread -lopenbabel
#LDFLAGS_USER =-L$(PAGMO_PATH)/src -lpagmo -lboost_serialization -lboost_system -lboost_thread -lpthread -lopenbabel
LDFLAGS_USER = -lpagmo -lboost_serialization -lboost_system -lboost_thread -lpthread -lopenbabel -lgsl -lgslcblas
$(BIN2): $(OBJS2)
	$(LINK_BINARY_CXX)

SRCS3 = cga_test.cpp sutton_chen_pot.cpp ga_utils.cpp
OBJS3 = ${SRCS3:.cpp=.o}

cga_test: $(OBJS3)
	g++ -o $@ $+ $(LDFLAGS) $(LDFLAGS_USER) -lgtest 

test: cga_test
