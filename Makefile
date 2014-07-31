SAMTOOLS_PATH = /gscratch/esci/dacb/samtools
SAMTOOLS_INCLUDE = $(SAMTOOLS_PATH)/include
SAMTOOLS_LIB = $(SAMTOOLS_PATH)/lib

CINCLUDE = -I. -I$(SAMTOOLS_INCLUDE)

CFLAGS	= $(COPTS) $(CINCLUDE) -g

# sources and object lists
SRCS	= bam2orf.c
OBJS	= $(SRCS:.c=.o)
PIC_OBJS= $(OBJS:.o=_PIC.o)
AR	= 
BINS	= bam2orf

default: all

all: $(BINS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

bam2orf: bam2orf.o $(AR)
	$(CC) $(CFLAGS) $< $(AR) -o $@ -L$(SAMTOOLS_LIB) -lbam -lz

Makefile.depend:
	touch Makefile.depend

depend: Makefile.depend
	makedepend -fMakefile.depend -- $(CFLAGS) -- $(SRCS)

clean: depend
	/bin/rm -rf $(BINS) $(OBJS) $(PIC_OBJS)

debug: $(BINS)
	gdb bam2orf -x gdb.run

include Makefile.depend
