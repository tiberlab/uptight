topdir := ../..

inputfile = parameterizer.in 

exec = $(topdir)/lib_uptight/PARAMETERIZER

.PHONY: run clean distclean


run: 
	@echo "Running test in $(notdir $(PWD)) ..."
	@echo ""
	@$(exec) < $(inputfile)


clean:
	@rm -rf *.dat


