all: bldall

bldall clean:
	cd  utils                &&  $(MAKE) $@
	cd  alg_lib              &&  $(MAKE) $@
	cd  tests                &&  $(MAKE) $@
