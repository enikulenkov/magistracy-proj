BUILD_HOME = ..
include $(BUILD_HOME)/Options

all: $(ALG_LIB)

CUDA_SRCS = matr_mult.cu
WRAP_SRCS = 

CUDA_OBJS = ${CUDA_SRCS:.cu=.o}

$(ALG_LIB): $(CUDA_OBJS)
	$(MAKE_LIB)
