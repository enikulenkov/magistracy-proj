BUILD_HOME=../../..
include $(BUILD_HOME)/Options

BIN1 = cublas_main
BINS= $(BIN1)

all: $(BINS)

SRCS= main.c
OBJS= ${SRCS:.c=.o}

LDFLAGS += -lcublas

$(BIN1): $(OBJS)
	$(LINK_BINARY_CUDA) $(BLAS_LIBS)
