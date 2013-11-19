BUILD_HOME=../..
include $(BUILD_HOME)/Options

BIN1 = sparse_matr_vec_mult
BINS= $(BIN1)

all: $(BINS)

SRCS= sparse_matr_vec_mult.c
OBJS= ${SRCS:.c=.o}

LDFLAGS += -lcusparse

$(BIN1): $(OBJS)
	$(LINK_BINARY)
