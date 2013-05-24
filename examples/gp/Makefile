CFLAGS = -std=gnu99 -O2
LFLAGS =
PKG_CONFIG_CFLAGS = 
PKG_CONFIG_LFLAGS = 
MODULE_OBJECTS = main_gp.o ga_gp.o programs.o

all : gp

gp : $(MODULE_OBJECTS)
	$(CC) $(CFLAGS) $(LFLAGS) $(MODULE_OBJECTS) -o gp -lm -lfgen $(PKG_CONFIG_LFLAGS)

main_gp.o : main_gp.c
	$(CC) -c $(CFLAGS) main_gp.c -o main_gp.o $(PKG_CONFIG_CFLAGS)

programs.o : programs.c
	$(CC) -c $(CFLAGS) programs.c -o programs.o $(PKG_CONFIG_CFLAGS)

ga_gp.o : ga_gp.c
	$(CC) -c $(CFLAGS) ga_gp.c -o ga_gp.o $(PKG_CONFIG_CFLAGS)

clean :
	rm -f $(MODULE_OBJECTS)
	rm -f gp

