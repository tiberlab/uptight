
hostos     = linux 

CXX      = g++


CXXFLAGS = -fPIC -I ../../../include/utils/UPTIGHT

# shared library extension
shared_libext := so
ifeq ($(findstring cygwin,$(hostos)),cygwin)
  shared_libext := dll
endif



all: libuptightc.$(shared_libext)


libuptightc.$(shared_libext):  UptWrapper.o
	@echo "Linking libuptightc.$(shared_libext) ..."
	@$(CXX) -shared -o $@ $< \
		-Wl,--version-script,export.map \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libuptight.a \
		-Wl,-Bdynamic

%.o: %.C
	@echo "Compiling $< ..."
	@$(CXX) $(CXXFLAGS) -c -o $@ $<


clean:
	@echo "cleaning up ./src/utils/UPTIGHT"
	@rm -f *.o

distclean: clean



