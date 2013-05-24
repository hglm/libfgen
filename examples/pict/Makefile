
CFLAGS = -std=gnu99 -O2
LFLAGS =
PKG_CONFIG_CFLAGS = `pkg-config --cflags gtk+-3.0`
PKG_CONFIG_LFLAGS = `pkg-config --libs gtk+-3.0`
MODULE_OBJECTS = main.o gtk.o ga.o block4x4.o decode.o

all : pict

pict : $(MODULE_OBJECTS)
	$(CC) $(CFLAGS) $(LFLAGS) $(MODULE_OBJECTS) -o pict -lm -lfgen $(PKG_CONFIG_LFLAGS)


main.o : main.c
	$(CC) -c $(CFLAGS) main.c -o main.o $(PKG_CONFIG_CFLAGS)

gtk.o : gtk.c
	$(CC) -c $(CFLAGS) gtk.c -o gtk.o $(PKG_CONFIG_CFLAGS)

ga.o : ga.c
	$(CC) -c $(CFLAGS) ga.c -o ga.o $(PKG_CONFIG_CFLAGS)

decode.o : decode.c
	$(CC) -c $(CFLAGS) decode.c -o decode.o

block4x4.o : block4x4.c
	$(CC) -c $(CFLAGS) block4x4.c -o block4x4.o $(PKG_CONFIG_CFLAGS)

clean :
	rm -f $(MODULE_OBJECTS)
	rm -f pict

dep:
	rm -f .depend
	make .depend

.depend:
	echo '# Module dependencies' >>.depend
	$(CC) -MM $(patsubst %.o,%.c,$(MODULE_OBJECTS)) >>.depend

include .depend

