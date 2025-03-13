# $Id: Makefile 564 2007-11-23 09:02:50Z maufder $

topdir := ../..

inputfile = dr_rel 

exec = $(topdir)/lib_uptight/UPT_TEST6

.PHONY: run clean distclean


run: 
	@echo "Running test in $(notdir $(PWD)) ..."
	@echo ""
	@$(exec) < $(inputfile)


clean:
	@rm -rf *.dat


