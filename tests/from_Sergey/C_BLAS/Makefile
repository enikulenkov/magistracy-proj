BUILD_HOME = ../../..
include $(BUILD_HOME)/Options

BINS = cblas_main

all: $(BINS)

SRCS = main.c

OBJS = ${SRCS:.c=.o}

cblas_main: $(OBJS)
	$(LINK_BINARY) $(BLAS_LIBS) 
