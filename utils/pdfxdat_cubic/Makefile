gCC=gcc
CFLAGS=-lm -g

SOURCES=pdfxdat_cubic.c
TARGETS=$(SOURCES:.c=)
HEADERS=
HSOURCES=$(HEADERS:.h=.c)

all : $(TARGETS)

.PHONY : all

$(TARGETS) : $(SOURCES) $(HEADERS) $(HSOURCES)
	$(CC)  -o $@ $(@).c  $(HSOURCES) $(CFLAGS)

clean:
	rm *.o
cleanall:
	rm  $(TARGETS)
