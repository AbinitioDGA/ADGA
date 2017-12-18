all: binR abinitiodga make_config 

make_config:
	cd make_configs/;./make_config_auto.sh

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

abinitiodga: make_config
	cd src/; make

clean: make_config
	cd src/; make clean

pristine: make_config
	cd src/; make pristine
	rmdir bin
