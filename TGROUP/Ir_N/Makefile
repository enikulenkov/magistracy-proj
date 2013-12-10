CFLAGS = -g -fdollar-ok
BIN = sc_qdm_n

all: clean $(BIN)

$(BIN): SC_QDM_N.f90
	gfortran -o $@ $(CFLAGS) $+

clean:
	rm -f $(BIN)
