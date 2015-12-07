# CVS info: $Id: Makefile,v 1.30 2008-12-14 22:06:58 mclamp Exp $

# For gcc
CC=cc

#CFLAGS= -g -Wall `pkg-config --cflags gtk+-2.0`
CFLAGS= -O3
#CFLAGS= -g
#CFLAGS=-g -O -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes `pkg-config --cflags gtk+-2.0`
#CFLAGS=-g -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes
LDFLAGS= -lm
#MYSQLCFLAGS=${shell mysql_config --cflags}
#MYSQLLIBS=${shell mysql_config --libs}
#INCDIR=/home/radon00/mclamp/cvs/ensc-core/src/
#ENSFLAGS =  -g -Wall -I$(INCDIR)DataModel -I$(INCDIR)DataAdaptors -I$(INCDIR)Util -I$(INCDIR)DataModel/Mapper -I$(INCDIR)DataModel/Basic -I$(INCDIR)Compara/DataAdaptors -I$(INCDIR)Compara/DataModel $(MYSQLCFLAGS)

all: random_corr test_fasta test_pwm test_pimer test_random_pos test_all gff_random_coverage test_pi gff_pi_coverage test_corr exon_pi_coverage test_fft gff_pimer_coverage gff_pimer_mercov gff_pimer_mercov_HMRD fasta2reads poghash scan_transfac parse_transout read_region gff_aligncov tfscan

clean:
	\rm -f *.o random_corr test_fasta test_pimer test_random_pos test_all gff_random_coverage test_pi gff_pi_coverage test_corr exon_pi_coverage test_fft gff_pimer_coverage gff_pimer_mercov gff_pimer_mercov_HMRD fasta2reads poghash scan_transfac parse_transout scan_transfac_gff read_region gff_aligncov tfscan test_pwm

##################################################

pog_utils.o: pog_utils.c pog_utils.h
	$(CC) $(CFLAGS) -c pog_utils.c

chromosome.o: chromosome.c
	$(CC) $(CFLAGS) -c chromosome.c

pimer.o: pimer.c
	$(CC) $(CFLAGS) -c pimer.c

fasta.o: fasta.c
	$(CC) $(CFLAGS) -c fasta.c

chain.o: chain.c
	$(CC) $(CFLAGS) -c chain.c

fasta2reads: fasta2reads.c fasta.o pog_utils.o
	$(CC) $(CFLAGS) -o fasta2reads fasta2reads.c fasta.o pog_utils.o ${LDFLAGS}

random_corr: random_corr.c fasta.o gff.o pimer.o chromosome.o pog_utils.o
	$(CC) $(CFLAGS) -o random_corr random_corr.c fasta.o gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

test_corr: test_corr.c fasta.o gff.o pimer.o chromosome.o pog_utils.o
	$(CC) $(CFLAGS) -o test_corr test_corr.c fasta.o gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

test_fasta: test_fasta.c fasta.o gff.o pimer.o chromosome.o pog_utils.o
	$(CC) $(CFLAGS) -o test_fasta test_fasta.c fasta.o gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

test_pwm: test_pwm.c fasta.o gff.o pimer.o chromosome.o pog_utils.o
	$(CC) $(CFLAGS) -o test_pwm test_pwm.c fasta.o gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

test_pimer: test_pimer.c gff.o pimer.o chromosome.o pog_utils.o
	$(CC) $(CFLAGS) -o test_pimer test_pimer.c gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

test_pi: test_pi.c gff.o pimer.o chromosome.o pog_utils.o
	$(CC) $(CFLAGS) -o test_pi test_pi.c gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

gff_random_coverage: gff_random_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o gff_random_coverage gff_random_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

exon_pi_coverage: exon_pi_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o exon_pi_coverage exon_pi_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

gff_pi_coverage: gff_pi_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o gff_pi_coverage gff_pi_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

gff_pimer_coverage: gff_pimer_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o gff_pimer_coverage gff_pimer_coverage.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

gff_pimer_mercov: gff_pimer_mercov.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o gff_pimer_mercov gff_pimer_mercov.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

gff_pimer_mercov_HMRD: gff_pimer_mercov_HMRD.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o gff_pimer_mercov_HMRD gff_pimer_mercov_HMRD.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

test_random_pos: test_random_pos.c gff.o pimer.o chromosome.o pog_utils.o 
	$(CC) $(CFLAGS) -o test_random_pos test_random_pos.c gff.o pimer.o chromosome.o pog_utils.o ${LDFLAGS}

test_fft: test_fft.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o test_fft test_fft.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

test_all: test_all.c gff.o pimer.o chromosome.o pog_utils.o fasta.o
	$(CC) $(CFLAGS) -o test_all test_all.c gff.o pimer.o chromosome.o pog_utils.o fasta.o ${LDFLAGS}

poghash: poghash.c pog_utils.o fasta.o 
	$(CC) $(CFLAGS) -o poghash poghash.c pog_utils.o fasta.o ${LDFLAGS}

scan_transfac: scan_transfac.c gff.o pog_utils.o fasta.o  chromosome.o
	$(CC) $(CFLAGS) -o scan_transfac scan_transfac.c gff.o pog_utils.o fasta.o chromosome.o ${LDFLAGS} 

tfscan: tfscan.c gff.o pog_utils.o fasta.o  chromosome.o
	$(CC) $(CFLAGS) -o tfscan tfscan.c gff.o pog_utils.o fasta.o chromosome.o ${LDFLAGS} 

parse_transout: parse_transout.c gff.o pog_utils.o fasta.o  chromosome.o
	$(CC) $(CFLAGS) -o parse_transout parse_transout.c gff.o pog_utils.o fasta.o chromosome.o ${LDFLAGS}

scan_transfac_gff: scan_transfac_gff.c gff.o pog_utils.o fasta.o  chromosome.o
	$(CC) $(CFLAGS) -o scan_transfac_gff scan_transfac_gff.c gff.o pog_utils.o fasta.o chromosome.o ${LDFLAGS}

gff_aligncov: gff_aligncov.c gff.o pog_utils.o fasta.o  chromosome.o
	$(CC) $(CFLAGS) -o gff_aligncov gff_aligncov.c gff.o pog_utils.o fasta.o chromosome.o ${LDFLAGS}

read_region: read_region.c gff.o pog_utils.o fasta.o  chromosome.o
	$(CC) $(CFLAGS) -o read_region read_region.c gff.o pog_utils.o fasta.o chromosome.o ${LDFLAGS}
##################################################
