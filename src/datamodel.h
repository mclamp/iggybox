#ifndef DATAMODEL_H
#define DATAMODEL_H

#define HUMAN_CHRNUM 23
#define DEBUG 1
#define NO_DEBUG 0 
#define ORGS 20


typedef struct _Thresh {
  char **ids;
  int  *thresh;
  int  *num;
  int count;
} Thresh;

typedef struct _Pwm {
  char *name;
  int   len;
  double *vals;
  double *inf;
} Pwm;

typedef struct _OrgHash {
  char **orgs;
  int  *vals1;
  int  *vals2;
  int  *vals3;
  int  mdcount;
  int  no_md;
  int  count;
  int *xcov;
  int  *cov;
  char **mercov;
} OrgHash;

typedef struct _Chromosome {
  char *name;
  unsigned long long int   length;
} Chromosome;

typedef struct _Genome {
	char *species;
	char *version;
	int   chrnum;
	Chromosome **chr;
	unsigned long long int   length;
} Genome;

typedef struct _GenomePosition {
	Chromosome *chr;
	unsigned long long int pos;
} GenomePosition;

typedef struct _Pimer {
  char  *chr;
  int    start;
  int    end;
  double logodds;
} Pimer;

typedef struct _pi {
	char *chr;
	int   start;
	int   end;
	float aval;
	float cval;
	float gval;
	float tval;
	float tree_length;
	double logodds;
} Pi;

typedef struct _Sequence {
  char    *name;
  char    *id;
  char    *sequence;
  int     length;
  int     start;
  int     end;
  int     selected;
  struct _Sequence *next;
} Sequence;

typedef struct _SeqArray {
  int len;
  Sequence **seq;
} SeqArray;

typedef struct _PimerArray {
  int len;
  Pimer **pimer;
} PimerArray;

typedef struct _Gff {
  char *seqname;
  char *source;
  char *feature;
  int start;
  int end;
  int strand;
  float score;
  int phase;
  char *hseqname;
  int hstart;
  int hend;
  int hstrand;
  struct _Gff *next;
} Gff;


typedef struct _Block {
  int position;
  int length;
  int count;
  struct _Block *next;
} Block;



Sequence    *make_Sequence(void);

void        *Sequence_set_name(Sequence *seq,char *name);
void        *Sequence_set_id(Sequence *seq,char   *id);
void        *Sequence_set_seqstring(Sequence *seq,char *str);
void        Sequence_add_Sequence(Sequence *seq,Sequence *newseq);
Sequence    *Sequence_get_Sequence_at(Sequence *list,int pos,int max);
char        Sequence_get_char_at(Sequence *seq,int position);
int         Sequence_count_Sequences(Sequence *list);
void        get_consensus_sequence(Sequence *, int *consensus[]);


Sequence *get_consensus_Sequence(Sequence *head);
int Sequence_max_sequence_length(Sequence *list);
Block **find_Blocks(Sequence *head);


extern char *substring(char *seq,int start, int end);

extern void  pog_free( void*ptr, char *desc, int debug);

extern Genome *fetch_human_genome();

extern GenomePosition *get_random_genome_position(Genome *genome, int offset);
extern void           free_random_genome_position(GenomePosition *pos);

extern Pimer          *read_pimer(FILE *file);
extern void            free_pimer(Pimer *pimer);

extern Pimer         **read_pimers_test(GenomePosition *pos,int len);
extern void            free_pimers_test(Pimer **pis, int len);

extern Pi             *read_pi(FILE *file);
extern Pi            **read_pis_test(GenomePosition *pos, int len);
extern void            free_pis_test(Pi **pis, int len);

extern Sequence       *read_fasta(FILE *file);
extern Sequence       *read_fasta_test(GenomePosition *pos,int len);
extern void            free_seqptr(Sequence *seq);

extern Gff            *read_gff(FILE *file);
extern void            print_gff(FILE *file,Gff *gff);
extern void            free_gffptr(Gff *gff);

extern Gff            *read_genome_position_gff(GenomePosition *pos, int len);

extern Chromosome     *find_chromosome_by_name(Genome *genome, char *name);


#endif /* DATAMODEL_H */

