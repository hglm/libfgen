
CC = gcc
CCPLUSPLUS = g++

# Uncomment the following 5 lines for normal compilation.
CFLAGS = -std=gnu99 -Wall -Ofast
CPPFLAGS = -Wall -Ofast
LIBFGEN = -lfgen
LIBFGENPP = -lfgenpp
LIBS = -lm
# Uncomment the following 5 lines for compilation with profiling/debugging.
#CFLAGS = -std=gnu99 -Wall -O -ggdb -p
#CPPFLAGS = -Wall -O -pg
#LIBFGEN = -lfgen_pg
#LIBFGENPP = -lfgenpp_pg
#LIBS = -lm -lpthread

LIB_DIR = ../. 
C_MODULES = *.c *.cpp

all : ones_bin pi_pso_bin func_bin tsp_bin tsp_archipelago_bin \
	  ffit_rosenbrock_bin ffit_rastrigin_bin rastrigin_double_bin tsp_archipelago_cpp_bin

ones_bin : ones.c
	$(CC) $(CFLAGS) -DFGEN_DLL ones.c -L$(LIB_DIR) -o ones_bin $(LIBFGEN) $(LIBS)

pi_pso_bin : pi_pso.c
	$(CC) $(CFLAGS) -DFGEN_DLL pi_pso.c -L$(LIB_DIR) -o pi_pso_bin $(LIBFGEN) $(LIBS)

func_bin : func.c
	$(CC) $(CFLAGS) -DFGEN_DLL func.c -L$(LIB_DIR) -o func_bin $(LIBFGEN) $(LIBS)

tsp_bin : tsp.c
	$(CC) $(CFLAGS) -DFGEN_DLL tsp.c -L$(LIB_DIR) -o tsp_bin $(LIBFGEN) $(LIBS)

tsp_archipelago_bin : tsp_archipelago.c
	$(CC) $(CFLAGS) -DFGEN_DLL tsp_archipelago.c -L$(LIB_DIR) -o tsp_archipelago_bin $(LIBFGEN) $(LIBS)

ffit_rosenbrock_bin : ffit_rosenbrock.c
	$(CC) $(CFLAGS) -DFGEN_DLL ffit_rosenbrock.c -L$(LIB_DIR) -o ffit_rosenbrock_bin $(LIBFGEN) $(LIBS)

ffit_rastrigin_bin : ffit_rastrigin.c
	$(CC) $(CFLAGS) -DFGEN_DLL ffit_rastrigin.c -L$(LIB_DIR) -o ffit_rastrigin_bin $(LIBFGEN) $(LIBS)

rastrigin_double_bin : rastrigin_double.c
	$(CC) $(CFLAGS) -DFGEN_DLL rastrigin_double.c -L$(LIB_DIR) -o rastrigin_double_bin $(LIBFGEN) $(LIBS)

tsp_archipelago_cpp_bin : tsp_archipelago_cpp.cpp
	$(CCPLUSPLUS) $(CPPFLAGS) -Wall -O2 -DFGEN_DLL tsp_archipelago_cpp.cpp -L$(LIB_DIR) -o tsp_archipelago_cpp_bin $(LIBFGENPP) $(LIBFGEN) $(LIBS)

clean :
	rm -f ones_bin
	rm -f pi_pso_bin
	rm -f func_bin
	rm -f tsp_bin
	rm -f tsp_archipelago_bin
	rm -f ffit_rosenbrock_bin
	rm -f ffit_rastrigin_bin
	rm -f rastrigin_double_bin
	rm -f tsp_archipelago_cpp_bin

dep:
	rm -f .depend
	make .depend

.depend:
	echo '# Module dependencies' >>.depend
	$(CC) -MM $(C_MODULES) >>.depend

include .depend

