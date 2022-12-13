HOMEPATH=$(PWD)
MODLOC=/home/zpiskulich_sta/privatemodules

wham_derivative:
	@echo "Making wham derivative code"
	cp modulefiles/wham-derivative.lua $(MODLOC)
	@echo "prepend_path('PATH', \"$(HOMEPATH)/bin\")" >> $(MODLOC)/wham_derivative.lua
	mkdir -p bin/
	touch bin/test
	rm bin/*
	ln -s $(HOMEPATH)/wham_derivatives/preprocess_wham.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/wham_derivatives/wham_class.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/wham_derivatives/wham_main.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/wham_derivatives/conv_distfiles.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/wham_derivatives/conv_metadata.py $(HOMEPATH)/bin/
	ln -s $(HOMEPATH)/wham_derivatives/pull_energies.py $(HOMEPATH)/bin/
	chmod 777 bin/*
