#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libconfig.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/sysctl.h>

#define BOOL char
#define TRUE 1
#define FALSE 0

#define PROBABILITY double

#define VERBOSE

void *ALLOC(size_t size);
  

// states and characters (emissions) are always considered to be numbered starting at 0

typedef struct s_r_s {
  int state_from; // which states you're allowed to be in at this position
  int state_to;
  struct s_r_s *next_range;
} state_range_struct;

typedef struct {
  int position_from; // which positions along the input sequence this corresponds to
  int position_to;
  state_range_struct *state_ranges;
} fixed_states_struct;

typedef struct {
  char **state_names; // n_states long array of pointers to strings containing state names
  int n_states;
  PROBABILITY *initial_probs;  // silent_states_begin array of probabilities of beginning in each non-silent state
  PROBABILITY *transition_matrix;  // n_states by n_states array of transition probabilities
                                   // consider row index to be "from" state, column index to be "to" state
//  PROBABILITY *transition_matrix_b;  // split this out to allow "fixed" states along a sequence (i.e. enforcement of a certain state at a certain position)
  PROBABILITY *emission_matrix;  // alphabet_length by n_states array of emission probabilities
                                 // consider row index to be state, column index to be alphabet index

  int silent_states_begin; // state index of first silent state.  setting equal to n_states means no silent states.  silent states are presumed to be ordered correctly for forward algorithm; backward algorithm will traverse them in the reverse order.

  char *alphabet;
  int alphabet_length;
//  bool case_sensitive; // is alphabet case sensitive?
  BOOL case_sensitive; // is alphabet case sensitive?

  int **parents;
  int *n_parents;
  int *first_silent_parent;

  int **children;
  int *n_children;
  int *first_silent_child;

  FILE *output;

  fixed_states_struct *fixed_states; // list of postions -> states restrictions for fixing input sequence positions to be in certain states
  int n_fixed_states; // how many fixed position -> state mappings we've got
  BOOL *fixed_state_positions; // vector indicating which positions have state restrictions
} model_def_struct;


typedef struct {
  char *seq; // array of output characters, represented by their index in alphabet in ASCII characters (i.e. letter 0 is chr(0))
  long len;
} sequence_struct;


typedef void(*a0k_func)(model_def_struct *);


typedef struct {
  model_def_struct *model_def;
  sequence_struct *sequence;
  PROBABILITY *table;
  PROBABILITY *s;
} thread_wrapper_struct;


// these functions will be helpful if I ever change how the matrices are constructed..
// from and to are given in state number
PROBABILITY fetch_transition_prob(model_def_struct *model_def, int from, int to);

void set_transition_prob(model_def_struct *model_def, int from, int to, PROBABILITY p);

// I don't think these need to be exposed for external usage
/*
PROBABILITY fetch_transition_prob_f(model_def_struct *model_def, int from, int to);

void set_transition_prob_f(model_def_struct *model_def, int from, int to, PROBABILITY p);

PROBABILITY fetch_transition_prob_b(model_def_struct *model_def, int from, int to);

void set_transition_prob_b(model_def_struct *model_def, int from, int to, PROBABILITY p);
*/

// state is given in state number, chr is given in character number
PROBABILITY fetch_emission_prob(model_def_struct *model_def, int state, int chr);

void set_emission_prob(model_def_struct *model_def, int state, int chr, PROBABILITY p);

void update_silent_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int row_index, BOOL forward);


void update_normal_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int row_index, BOOL forward);

PROBABILITY normalize_row(PROBABILITY *row, int n_states, int total_states);

/* INPUTS:
   model_def: struct containing definition of the model
   sequence: struct containing sequence to train the model on
   OUTPUTS:
   table: forward table
          this has to be filled out row-wise, because of paging issues.  hence, row index is sequence position and column index is state
   s: array of s_i probability scaling factors
*/
void forward(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, PROBABILITY *s);


/* INPUTS:
   model_def: struct containing definition of the model
   sequence: struct containing sequence to train the model on
   s: array of s_i probability scaling factors
   OUTPUTS:
   table: backward table
          this has to be filled out row-wise, because of paging issues.  hence, row index is sequence position and column index is state
          this table is stored in reverse (i.e. row 0 contains the last column of the backwards table), also for paging reasons
*/
void backward(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *s, PROBABILITY *table);


PROBABILITY posterior_decoding(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *sb, PROBABILITY *sr, int obs_index, int state);


void calc_sr(PROBABILITY *sf, PROBABILITY *sb, int len, PROBABILITY *sr);


void print_forward_table(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, PROBABILITY *s, int n);


void print_backward_table(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, PROBABILITY *s, int n);


void print_posterior(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *s, int states_to_include, int n);


void print_most_probable_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *s);


PROBABILITY A_kl(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, int n_seqs, int k, int l);


PROBABILITY E_kb(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, PROBABILITY **s, int n_seqs, int k, char b);

void print_transition_matrix_as_dot(model_def_struct *model_def);


void print_transition_matrix(model_def_struct *model_def);


PROBABILITY log_likelihood(PROBABILITY **s, sequence_struct **sequence, int n_seqs);


void baum_welch(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, PROBABILITY **s, int n_seqs, PROBABILITY delta, a0k_func update_a0k);


model_def_struct *initialize_model(char *filename, fixed_states_struct* fixed_states, int n_fixed_states);

void find_parents_and_children(model_def_struct *model_def);

int read_sequence(char *filename, sequence_struct ***sequence_ptr);


void print_initial_probs(model_def_struct *model_def);


BOOL verify_model(model_def_struct *model_def);


BOOL anything_to_update(model_def_struct *model_def);

int find_num_cpus();


typedef struct {
  int index;
  PROBABILITY value;
} sortable_pair;
