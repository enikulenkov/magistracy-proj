BUILD_HOME = ..
include $(BUILD_HOME)/Options

all: $(UTILS_LIB)

SRCS = mat_file_parser.c utils.c ts_util.c timer.c

OBJS = ${SRCS:.c=.o}

$(UTILS_LIB): $(OBJS)
	$(MAKE_LIB)
