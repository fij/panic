INCFLAGS=-I.
LINKFLAGS=-lm
CFLAGS=-O3 $(INCFLAGS) $(LINKFLAGS)
W_CFLAGS=-O3 -Wall $(INCFLAGS) $(LINKFLAGS)

all: demo

dc: demo crunch
	
clean:
	rm -f core *~ 

d: demo
c: crunch

demo: 
	gcc sd.c $(CFLAGS) -L/usr/X11R6/lib -lX11 -o sd

crunch:
	gcc sd_crunch.c $(CFLAGS) -o sd_crunch

#----------------------------------

dw:
	gcc sd.c $(W_CFLAGS) -L/usr/X11R6/lib -lX11 -o sd

cw: 
	gcc sd_crunch.c $(W_CFLAGS) -o sd_crunch


# DO NOT DELETE THIS LINE


