#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "sam.h"

#define MAX_LINE 102400

#define CDS 1
#define tRNA 2
#define rRNA 3

typedef struct gene {
	char *locus;
	uint8_t type;
	int start;
	int end;
	char *locus_tag;
	long reads;
} gene_t;

typedef struct gff {
	char *filename;
	int genes;
	gene_t *gene_features;
} gff_t;

char *sndup(const char *fmt, ...) {
	char    *buf;
	int     size = strlen(fmt) * sizeof(char);
	int     ret;
	va_list va_l;

	buf = (char *)malloc(size);

	while (1) {
		va_start(va_l, fmt);
		ret = vsnprintf(buf, size, fmt, va_l);
		va_end(va_l);

		if (ret > -1 && ret < size)
			return buf;

		if (ret > -1)
			size = ret + 1;
		else
			size *= 2;

		buf = (char *)realloc(buf, size);
	}
}

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);
	return 0;
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	gene_t *gf = (gene_t *)data;
	if ((int)pos >= gf->start && (int)pos <= gf->end)
		gf->reads += n;
	return 0;
}

gff_t *gffopenreadclose(char *filename) {
	FILE *fp;
	gff_t *gf;
	char line[MAX_LINE], *last, *locus, *type, *start, *end, *attrs, *attr, *locus_tag;
	gene_t *gene;
	int lines;

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "unable to open GFF file '%s'\n", filename);
		return NULL;
	}

	gf = (gff_t *)malloc(sizeof(gff_t));
	gf->filename = strdup(filename);
	gf->genes = 0;
	gf->gene_features = NULL;

	while (fgets(line, MAX_LINE, fp) != NULL) {
		++lines;
		if (line[0] == '#') continue;

		locus = strtok_r(line, "\t", &last);
		strtok_r((char *)NULL, "\t", &last);		// ignore field
		type = strtok_r((char *)NULL, "\t", &last);
		if (strcmp(type, "gene") == 0) {
			start = strtok_r((char *)NULL, "\t", &last);
			end = strtok_r((char *)NULL, "\t", &last);
			strtok_r((char *)NULL, "\t", &last);		// ignore field
			strtok_r((char *)NULL, "\t", &last);		// ignore field
			strtok_r((char *)NULL, "\t", &last);		// ignore field
			attrs = strtok_r((char *)NULL, "\t\n", &last);

			gf->gene_features = (gene_t *)realloc(gf->gene_features, (gf->genes + 1) * sizeof(gene_t));
			gene = &gf->gene_features[gf->genes];

			gene->locus = strdup(locus);
			
			gene->start = atoi(start);
			gene->end = atoi(end);
			gene->reads = 0;
			gene->type = 0;

			// parse attribute line for ID
			attr = strtok_r(attrs, ";", &last);
			while(attr != NULL) {
				if (strncmp(attr, "ID=", 3) == 0) {
					gene->locus_tag = sndup("%s", &attr[3]);
					break;
				}
				attr = strtok_r((char *)NULL, ";", &last);
			}
			if (attr == NULL) {
				gene->locus_tag = sndup("%s:%d-%d", gene->locus, gene->start, gene->end);
			}

			gf->genes++;
		} else if (strcmp(type, "CDS") == 0) {
			start = strtok_r((char *)NULL, "\t", &last);
			end = strtok_r((char *)NULL, "\t", &last);

			if (gene && gene->start == atoi(start) && gene->end == atoi(end))
				gene->type = CDS;
		} else
			continue;
	}

	fclose(fp);

	return gf;
}

int main(int argc, char *argv[]) {
	int g, read_length;
	gff_t *gf;
	samfile_t *bf;
	bam_index_t *bif;
	double reads_mapped_sum = 0;

	if (argc != 4) {
		fprintf(stderr, "Usage: %s <genome GFF> <sorted / indexed bam> <read length>\n", argv[0]);
		return 1;
	}

	gf = gffopenreadclose(argv[1]);
	if (!gf) {
		fprintf(stderr, "%s: an error occurred reading the GFF file '%s'\n", argv[0], argv[1]);
		return 1;
	}
	fprintf(stderr, "%s: found %d gene features in GFF file '%s'\n", argv[0], gf->genes, argv[1]);

	bf = samopen(argv[2], "rb", 0);
	if (!bf) {
		fprintf(stderr, "%s: an error occurred opening the BAM file '%s'\n", argv[0], argv[2]);
		return 1;
	}

	bif = bam_index_load(argv[2]);
	if (!bif) {
		fprintf(stderr, "%s: an error occurred opening the BAM index for BAM file '%s'\n", argv[0], argv[2]);
		return 1;
	}

	read_length = atoi(argv[3]);

	for (g = 0; g < gf->genes; ++g) {
		gene_t *gene = &gf->gene_features[g];
		char *region = sndup("%s:%d-%d", gene->locus, (gene->start < 0 ? 0 : gene->start), gene->end);
		int ref, start, end;
		bam_plbuf_t *buf;

		bam_parse_region(bf->header, region, &ref, &start, &end);
		if (ref < 0) {
			fprintf(stderr, "%s: invalid region '%s'\n", argv[0], region);
			return 1;
		}

		buf = bam_plbuf_init(pileup_func, gene); // initialize pileup
		// fix for pileup capping # of alignments to 8000
		bam_plp_set_maxcnt(buf->iter, 1000000);
		bam_fetch(bf->x.bam, bif, ref, start, end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_plbuf_destroy(buf);

		free(region);
	}
	bam_index_destroy(bif);
	samclose(bf);

	for (g = 0; g < gf->genes; ++g) {
		gene_t *gene = &gf->gene_features[g];
		if (gene->type == CDS)
			reads_mapped_sum += (double)gene->reads / (double)read_length;
	}
	reads_mapped_sum /= 1000000.;
	fprintf(stderr, "%s: millions of reads mapped to CDS = %.2f\n", argv[0], reads_mapped_sum);

	printf("#locus_tag\tRPKM\treads_mapped\n");
	for (g = 0; g < gf->genes; ++g) {
		gene_t *gene = &gf->gene_features[g];
		double rpkm = (double)gene->reads / (double)read_length / ((double)(gene->end - gene->start + 1) / 1000.) / reads_mapped_sum;
		printf("%s\t%f\t%f\n", gene->locus_tag, rpkm, (double)gene->reads / (double)read_length);
//		printf("%s\t%f\t%f\t%f\t%f\t%f\n", gene->locus_tag, rpkm, (double)gene->reads / (double)read_length, read_length, (gene->end -  gene->start + 1.) /1000., reads_mapped_sum);
	}

	return 0;
}
