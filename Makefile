#
# "make" compiles the library.
# "sudo make install" installs the library.
# "make profiling" makes the static library with profiling and debugging enabled.

CC = gcc
CPLUSPLUS = g++
SHARED_CFLAGS = -std=gnu99 -Wall -Ofast -fPIC -fvisibility=hidden -I. -DFGEN_SHARED -DFGEN_SHARED_EXPORTS
SHARED_CPLUSPLUSFLAGS = -Wall -Ofast -fPIC -fvisibility=hidden -I. -DFGEN_SHARED -DFGEN_SHARED_EXPORTS
PROFILING_CFLAGS = -std=gnu99 -Wall -ggdb -p -I.
PROFILING_CPLUSPLUSFLAGS = -Wall -ggdb -p -I.
SHARED_LIB_DIR = /usr/lib
INCLUDE_DIR = /usr/include
VERSION = 0.2.2
VERSION_MAJOR = 0

MODULES = bitstring.o parameters.o population.o \
	  random.o selection.o crossover.o mutation.o error.o \
	  gray.o ga.o decode.o seed.o ffit.o cache.o migration.o pso.o steady_state.o
FGENPP_MODULES = fgenpp.o
MODULES_PG = bitstring.p parameters.p population.p \
	  random.p selection.p crossover.p mutation.p error.p \
	  gray.p ga.p decode.p seed.p ffit.p cache.p migration.p pso.p steady_state.p
FGENPP_MODULES_PG = fgenpp.p

all : libfgen.so.$(VERSION) libfgenpp.so.$(VERSION)

libfgen.so.$(VERSION) : $(MODULES)
	$(CC) -shared -Wl,-soname,libfgen.so.$(VERSION_MAJOR) -fPIC -o libfgen.so.$(VERSION) $(MODULES) -lc -lm -lpthread

libfgenpp.so.$(VERSION) : $(FGENPP_MODULES)
	$(CPLUSPLUS) -shared -Wl,-soname,libfgenpp.so.$(VERSION_MAJOR) -fPIC -o libfgenpp.so.$(VERSION) $(FGENPP_MODULES)

install : $(SHARED_LIB_DIR)/libfgen.so.$(VERSION) $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION) $(INCLUDE_DIR)/fgen.h $(INCLUDE_DIR)/fgenpp.h

uninstall :
	rm -f $(SHARED_LIB_DIR)/libfgen.so.$(VERSION)
	rm -f $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION)
	rm -f $(SHARED_LIB_DIR)/libfgen.so
	rm -f $(SHARED_LIB_DIR)/libfgen.so.$(VERSION_MAJOR)
	rm -f $(SHARED_LIB_DIR)/libfgenpp.so
	rm -f $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION_MAJOR)
	rm -f $(INCLUDE_DIR)/fgen.h
	rm -f $(INCLUDE_DIR)/fgenpp.h

$(SHARED_LIB_DIR)/libfgen.so.$(VERSION) : libfgen.so.$(VERSION)
	install -m 0644 libfgen.so.$(VERSION) $(SHARED_LIB_DIR)/libfgen.so.$(VERSION)
	ln -sf $(SHARED_LIB_DIR)/libfgen.so.$(VERSION) $(SHARED_LIB_DIR)/libfgen.so
	ln -sf $(SHARED_LIB_DIR)/libfgen.so.$(VERSION) $(SHARED_LIB_DIR)/libfgen.so.$(VERSION_MAJOR)

$(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION) : libfgenpp.so.$(VERSION)
	install -m 0644 libfgenpp.so.$(VERSION) $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION)
	ln -sf $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION) $(SHARED_LIB_DIR)/libfgenpp.so
	ln -sf $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION) $(SHARED_LIB_DIR)/libfgenpp.so.$(VERSION_MAJOR)

$(INCLUDE_DIR)/fgen.h : fgen.h
	install -m 0644 fgen.h $(INCLUDE_DIR)/fgen.h

$(INCLUDE_DIR)/fgenpp.h : fgenpp.h
	install -m 0644 fgenpp.h $(INCLUDE_DIR)/fgenpp.h

.c.o :
	$(CC) -c $(SHARED_CFLAGS) $< -o $@

.cpp.o :
	$(CPLUSPLUS) -c $(SHARED_CPLUSPLUSFLAGS) $< -o $@

profiling : libfgen_pg.a libfgenpp_pg.a

libfgen_pg.a : $(MODULES_PG)
	ar r libfgen_pg.a $(MODULES_PG)

libfgenpp_pg.a : $(FGENPP_MODULES_PG)
	ar r libfgenpp_pg.a $(FGENPP_MODULES_PG)

.c.p :
	$(CC) -c $(PROFILING_CFLAGS) $< -o $@

.cpp.p :
	$(CPLUSPLUS) -c $(PROFILING_CPLUSPLUSFLAGS) $< -o $@

clean :
	rm -f $(MODULES)
	rm -f $(FGENPP_MODULES)
	rm -f libfgen.so.$(VERSION)
	rm -f libfgenpp.so.$(VERSION)
	rm -f $(MODULES_PG)
	rm -f $(FGENPP_MODULES_PG)
	rm -f libfgen_pg.a
	rm -f libfgenpp_pg.a

dep:
	rm -f .depend
	make .depend

.depend:
	echo '# Module dependencies' >>.depend
	$(CC) -MM $(patsubst %.o,%.c,$(MODULES)) >>.depend
	$(CC) -MM $(patsubst %.o,%.cpp,$(FGENPP_MODULES)) >>.depend
	for x in $(patsubst %.p,%.c,$(MODULES_PG)); do \
	$(CC) -MM -MT `echo $$x | sed s/\\\.c/\.p/` $$x >>.depend; \
	done; \
	for x in $(patsubst %.p,%.cpp,$(FGENPP_MODULES_PG)); do \
	$(CC) -MM -MT `echo $$x | sed s/\\\.cpp/\.p/` $$x >>.depend; \
	done

include .depend

