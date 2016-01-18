#include "bc.h"


void *ALLOC(size_t size) {
  void *foo;
  if(!(foo = malloc(size))) {
    fprintf(stderr, "malloc failed.\n");
    exit(0);
  }

  return foo;
}
  
// these functions will be helpful if I ever change how the matrices are constructed..
// from and to are given in state number
PROBABILITY fetch_transition_prob(model_def_struct *model_def, int from, int to) {
  return model_def->transition_matrix[from * model_def->n_states + to];
}

void set_transition_prob(model_def_struct *model_def, int from, int to, PROBABILITY p) {
  model_def->transition_matrix[from * model_def->n_states + to] = p;
}

// state is given in state number, chr is given in character number
PROBABILITY fetch_emission_prob(model_def_struct *model_def, int state, int chr) {
  // the emission matrix index is state number * number of characters in the alphabet + offset of the character observed here
  return model_def->emission_matrix[state * model_def->alphabet_length + chr];
}

void set_emission_prob(model_def_struct *model_def, int state, int chr, PROBABILITY p) {
  model_def->emission_matrix[state * model_def->alphabet_length + chr] = p;
}

void update_silent_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int row_index, BOOL forward) {
  PROBABILITY *row = table + (unsigned long)model_def->n_states * (unsigned long)row_index;
  int i, j;

  if (forward) {
    for (i = model_def->silent_states_begin; i < model_def->n_states; i++) {
      PROBABILITY sum = 0;
      for (j = 0; j < model_def->first_silent_parent[i]; j++) {
        int state = model_def->parents[i][j];
        PROBABILITY transition_prob = fetch_transition_prob(model_def, state, i);
        sum += row[state] * transition_prob;
      }
      row[i] = sum;
    }

    for (i = model_def->silent_states_begin + 1; i < model_def->n_states; i++) {
      PROBABILITY sum = 0;
      for (j = model_def->first_silent_parent[i]; j < model_def->n_parents[i]; j++) {
        int state = model_def->parents[i][j];
//        if (state >= i) break;  // I don't think this can happen, due to required ordering of silent states
        PROBABILITY transition_prob = fetch_transition_prob(model_def, state, i);
        sum += row[state] * transition_prob;
      }
      row[i] += sum;
    }

  } else {

    for (i = model_def->silent_states_begin; i < model_def->n_states; i++) {
      PROBABILITY sum = 0;
      for (j = 0; j < model_def->first_silent_child[i]; j++) {
        int state = model_def->children[i][j];
        PROBABILITY transition_prob = fetch_transition_prob(model_def, i, state);
        PROBABILITY emission_prob = fetch_emission_prob(model_def, state, sequence->seq[sequence->len - row_index - 1]);
        sum += row[state] * transition_prob * emission_prob;
      }
      row[i] = sum;
    }

    for (i = model_def->n_states - 2; i > model_def->silent_states_begin - 1; i--) {
      PROBABILITY sum = 0;
      for (j = model_def->first_silent_child[i]; j < model_def->n_children[i]; j++) {
        int state = model_def->children[i][j];
        PROBABILITY transition_prob = fetch_transition_prob(model_def, i, state);
        sum += row[state] * transition_prob;
      }
      row[i] += sum;
    }
  }
}

void update_normal_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int row_index, BOOL forward) {
  PROBABILITY *row = table + (unsigned long)model_def->n_states * (unsigned long)row_index;
  PROBABILITY *prev_row = table + (unsigned long)model_def->n_states * (unsigned long)(row_index - 1);
  PROBABILITY *next_row = table + (unsigned long)model_def->n_states * (unsigned long)(row_index + 1);
  int i, j;

  if (forward) {
    for (i = 0; i < model_def->silent_states_begin; i++) {
      PROBABILITY sum = 0;
      for (j = 0; j < model_def->n_parents[i]; j++) {
        int state = model_def->parents[i][j];
        PROBABILITY transition_prob = fetch_transition_prob(model_def, state, i);
        sum += prev_row[state] * transition_prob;
      }
      row[i] = fetch_emission_prob(model_def, i, sequence->seq[row_index]) * sum;
    }
  } else {
    for (i = 0; i < model_def->silent_states_begin; i++) {
      PROBABILITY sum = 0;
      for (j = 0; j < model_def->n_children[i]; j++) {
        int state = model_def->children[i][j];
        PROBABILITY transition_prob = fetch_transition_prob(model_def, i, state);
        PROBABILITY emission_prob = fetch_emission_prob(model_def, state, sequence->seq[sequence->len - row_index]);
        sum += prev_row[state] * transition_prob * emission_prob;
      }
      row[i] = sum;
    }
  }
}

PROBABILITY normalize_row(PROBABILITY *row, int n_states, int total_states) {
  PROBABILITY s = 0;
  int i;

//  return 1;

  for (i = 0; i < n_states; i++) {
    s += row[i];
  }
  for (i = 0; i < total_states; i++) {
    row[i] /= s;
  }

  return s;
}


void find_all_silent_children(model_def_struct *model_def, BOOL *buffer, int state, int *n_parents) {
  int i;

  if (model_def->first_silent_parent[state] >= model_def->n_parents[state]) { return; }

  for (i = model_def->first_silent_parent[state]; i < model_def->n_parents[state]; i++) {
    // make sure it's not already in the list, and if not, add it
    int parent = model_def->parents[state][i];

    if (!buffer[state * model_def->n_states + parent]) {
      buffer[state * model_def->n_states + parent] = TRUE;
      n_parents[state]++;
    }

    // it's silent; must check its parents, too
    find_all_silent_children(model_def, buffer, parent, n_parents);
  }
}


void find_all_silent_parents(model_def_struct *model_def, BOOL *buffer, int state, int *n_children) {
  int i, l;

  if (model_def->first_silent_parent[state] >= model_def->n_parents[state]) { return; }

  for (i = model_def->first_silent_parent[state]; i < model_def->n_parents[state]; i++) {
    // make sure it's not already in the list, and if not, add it
    int parent = model_def->parents[state][i];

    if (!buffer[parent * model_def->n_states + state]) {
      buffer[parent * model_def->n_states + state] = TRUE;
      n_children[parent]++;
    }

    for (l = 0; l < model_def->n_parents[parent]; l++) {
      if (!buffer[model_def->parents[parent][l] * model_def->n_states + parent]) {
        buffer[model_def->parents[parent][l] * model_def->n_states + parent] = TRUE;
        n_children[model_def->parents[parent][l]]++;
      }
    }

    // it's silent; must check its parents, too
    find_all_silent_parents(model_def, buffer, parent, n_children);
  }
}


/* INPUTS:
   model_def: struct containing definition of the model
   sequence: struct containing sequence to train the model on
   OUTPUTS:
   table: forward table
          this has to be filled out row-wise, because of paging issues.  hence, row index is sequence position and column index is state
   s: array of s_i probability scaling factors
*/
void forward(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, PROBABILITY *s) {

  // initialize first row
  // normal states need to be handled specially
  int i, j, k, l;
  for (i = 0; i < model_def->silent_states_begin; i++) {
    table[i] = model_def->initial_probs[i] * fetch_emission_prob(model_def, i, sequence->seq[0]);
  }

  // silent states can use the normal machinery
  update_silent_states(model_def, sequence, table, 0, TRUE);

  // the first weight is just the sum of the first column (row, in this implementation).  calculate, and normalize.
  s[0] = normalize_row(table, model_def->silent_states_begin, model_def->n_states);

  // first row initialization complete


  // fill out the rest of the table now
  for (i = 1; i < sequence->len; i++) {
    #ifdef VERBOSE
    if (i % 500 == 0) fprintf(stderr, "forward row %d\n", i);
    #endif

    // handle fixed state specifications entirely through parents and children lists
    int *n_parents, *n_parents_fixed, *first_silent_parent, *first_silent_parent_fixed;
    int **parents, **parents_fixed;
    BOOL *parents_matrix;
    BOOL handle_fixed_states = (model_def->n_fixed_states > 0) && model_def->fixed_state_positions[i];
    if (handle_fixed_states) {
      // there are fixed position ranges given that cover this position
      n_parents = model_def->n_parents;
      parents = model_def->parents;
      first_silent_parent = model_def->first_silent_parent;
      parents_matrix = ALLOC(sizeof(BOOL) * model_def->n_states * model_def->n_states);
      memset(parents_matrix, FALSE, sizeof(BOOL) * model_def->n_states * model_def->n_states);

      // build a new parents list to reflect the restrictions of the pinned positions
      parents_fixed = ALLOC(sizeof(int*) * model_def->n_states);
      n_parents_fixed = ALLOC(sizeof(int) * model_def->n_states);
      first_silent_parent_fixed = ALLOC(sizeof(int) * model_def->n_states);
      memset(parents_fixed, 0, sizeof(int*) * model_def->n_states);
      memset(n_parents_fixed, 0, sizeof(int) * model_def->n_states);
      memset(first_silent_parent_fixed, 0, sizeof(int) * model_def->n_states);

      for (j = 0; j < model_def->n_fixed_states; j++) {
        if ((i >= model_def->fixed_states[j].position_from) && (i <= model_def->fixed_states[j].position_to)) {
          // position i is covered in range j
          state_range_struct *state_range;
          for (state_range = model_def->fixed_states[j].state_ranges; state_range; state_range = state_range->next_range) {
            for (k = state_range->state_from; k <= state_range->state_to; k++) {
              // walk through all the states in the range and build a new parents list that includes only things from this range
              n_parents_fixed[k] = n_parents[k];
              parents_fixed[k] = parents[k];
              first_silent_parent_fixed[k] = model_def->first_silent_parent[k];
            }
          }

          for (k = model_def->silent_states_begin; k < model_def->n_states; k++) {
              n_parents_fixed[k] = n_parents[k];
              parents_fixed[k] = parents[k];
              first_silent_parent_fixed[k] = model_def->first_silent_parent[k];
          }
        }
      }

      model_def->n_parents = n_parents_fixed;
      model_def->parents = parents_fixed;
      model_def->first_silent_parent = first_silent_parent_fixed;
    }

    update_normal_states(model_def, sequence, table, i, TRUE);
    update_silent_states(model_def, sequence, table, i, TRUE);
    s[i] = normalize_row(table + (unsigned long)model_def->n_states * (unsigned long)i, model_def->silent_states_begin, model_def->n_states);

    if (handle_fixed_states) {
      model_def->n_parents = n_parents;
      model_def->parents = parents;
      model_def->first_silent_parent = first_silent_parent;
      free(n_parents_fixed);
      free(parents_fixed);
      free(first_silent_parent_fixed);
    }
  }
}


/* INPUTS:
   model_def: struct containing definition of the model
   sequence: struct containing sequence to train the model on
   s: array of s_i probability scaling factors
   OUTPUTS:
   table: backward table
          this has to be filled out row-wise, because of paging issues.  hence, row index is sequence position and column index is state
          this table is stored in reverse (i.e. row 0 contains the last column of the backwards table), also for paging reasons
*/
void backward(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *s, PROBABILITY *table) {

  // initialize first row
  // normal states need to be handled specially
  int i, j, k, l;

  s[sequence->len - 1] = model_def->silent_states_begin;
  for (i = 0; i < model_def->silent_states_begin; i++) {
    table[i] = 1.0 / s[sequence->len - 1];
  }

  // silent states can use the normal machinery
  update_silent_states(model_def, sequence, table, 0, FALSE);

  // first row initialization complete


  // fill out the rest of the table now
  for (i = 1; i < sequence->len; i++) {
    #ifdef VERBOSE
    if (i % 500 == 0) fprintf(stderr, "backward row %d\n", i);
    #endif
    int seq_pos = sequence->len - i - 1;

    // handle fixed state specifications entirely through parents and children lists
    int *n_children, *n_children_fixed, *first_silent_child, *first_silent_child_fixed;
    int **children, **children_fixed;
    BOOL *children_matrix;
    BOOL handle_fixed_states = (model_def->n_fixed_states > 0) && (model_def->fixed_state_positions[seq_pos] || model_def->fixed_state_positions[seq_pos + 1]);
    if (handle_fixed_states) {
      // there are fixed position ranges given that cover this position
      n_children = model_def->n_children;
      children = model_def->children;
      first_silent_child = model_def->first_silent_child;
      children_matrix = ALLOC(sizeof(BOOL) * model_def->n_states * model_def->n_states);
      memset(children_matrix, FALSE, sizeof(BOOL) * model_def->n_states * model_def->n_states);

      // build a new children list to reflect the restrictions of the pinned positions
      children_fixed = ALLOC(sizeof(int*) * model_def->n_states);
      n_children_fixed = ALLOC(sizeof(int) * model_def->n_states);
      first_silent_child_fixed = ALLOC(sizeof(int) * model_def->n_states);
      memset(children_fixed, 0, sizeof(int*) * model_def->n_states);
      memset(n_children_fixed, 0, sizeof(int) * model_def->n_states);
      memset(first_silent_child_fixed, 0, sizeof(int) * model_def->n_states);

      for (j = 0; j < model_def->n_fixed_states; j++) {
        if ((seq_pos + 1 >= model_def->fixed_states[j].position_from) && (seq_pos + 1 <= model_def->fixed_states[j].position_to)) {
          // position seq_pos + 1 is covered in range j
          state_range_struct *state_range;
          for (state_range = model_def->fixed_states[j].state_ranges; state_range; state_range = state_range->next_range) {
            for (k = state_range->state_from; k <= state_range->state_to; k++) {
              // walk through all the states in the range and build a new children list that includes only things from this range
              for (l = 0; l < model_def->first_silent_parent[k]; l++) {
                if (!children_matrix[model_def->parents[k][l] * model_def->n_states + k]) {
                  children_matrix[model_def->parents[k][l] * model_def->n_states + k] = TRUE;
                  n_children_fixed[model_def->parents[k][l]]++;
                }
              }
            }
          }

          for (k = model_def->silent_states_begin; k < model_def->n_states; k++) {
            // copy over the existing silent state dependencies
            n_children_fixed[k] = n_children[k];

            find_all_silent_parents(model_def, children_matrix, k, n_children_fixed);

            for (l = 0; l < n_children[k]; l++) {
              if (!children_matrix[k * model_def->n_states + model_def->children[k][l]]) {
                children_matrix[k * model_def->n_states + model_def->children[k][l]] = TRUE;
              }
            }
          }
        }

        if ((seq_pos >= model_def->fixed_states[j].position_from) && (seq_pos <= model_def->fixed_states[j].position_to)) {
          // if I'm at the position, I have to restrict the paths of the silent states
          for (k = 0; k < model_def->silent_states_begin; k++) {
            // copy over the existing normal state dependencies
            n_children_fixed[k] = n_children[k];

            for (l = 0; l < n_children[k]; l++) {
              if (!children_matrix[k * model_def->n_states + model_def->children[k][l]]) {
                children_matrix[k * model_def->n_states + model_def->children[k][l]] = TRUE;
              }
            }
          }

          for (k = model_def->silent_states_begin; k < model_def->n_states; k++) {
            n_children_fixed[k] = 0;

            // reset the existing silent state dependencies
            for (l = 0; l < n_children[k]; l++) {
              if (children_matrix[k * model_def->n_states + model_def->children[k][l]]) {
                children_matrix[k * model_def->n_states + model_def->children[k][l]] = FALSE;
              }
            }
          }

          state_range_struct *state_range;
          for (state_range = model_def->fixed_states[j].state_ranges; state_range; state_range = state_range->next_range) {
            for (k = state_range->state_from; k <= state_range->state_to; k++) {
              // walk up the chain of silent states to figure out the children list
              find_all_silent_parents(model_def, children_matrix, k, n_children_fixed);
            }
          }

        }
      }

      int counter = 0;
      for (j = 0; j < model_def->n_states; j++) {
        BOOL found_first_silent_child = FALSE;
        children_fixed[j] = ALLOC(sizeof(int) * n_children_fixed[j]);
        memset(children_fixed[j], 0, sizeof(int) * n_children_fixed[j]);
        for (k = 0; k < model_def->n_states; k++) {
          if (children_matrix[j * model_def->n_states + k]) {
            if (!found_first_silent_child && (k >= model_def->silent_states_begin)) {
              first_silent_child_fixed[j] = counter;
              found_first_silent_child = TRUE;
            }
            children_fixed[j][counter] = k;
            counter++;
          }
        }
        if (!found_first_silent_child) {
          first_silent_child_fixed[j] = counter;
        }
        counter = 0;
      }

      model_def->n_children = n_children_fixed;
      model_def->children = children_fixed;
      model_def->first_silent_child = first_silent_child_fixed;
    }

    update_normal_states(model_def, sequence, table, i, FALSE);
    update_silent_states(model_def, sequence, table, i, FALSE);
    s[seq_pos] = normalize_row(table + (unsigned long)model_def->n_states * (unsigned long)i, model_def->silent_states_begin, model_def->n_states);

    if (handle_fixed_states) {
      model_def->n_children = n_children;
      model_def->children = children;
      model_def->first_silent_child = first_silent_child;
      free(n_children_fixed);
      free(children_fixed);
      free(children_matrix);
      free(first_silent_child_fixed);
    }
  }
}


PROBABILITY fetch_forward_prob(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int obs_index, int state) {
  return table[(unsigned long)model_def->n_states * (unsigned long)obs_index + state];
}


PROBABILITY fetch_backward_prob(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int obs_index, int state) {
  return table[(unsigned long)model_def->n_states * (unsigned long)(sequence->len - obs_index - 1) + state];
}


PROBABILITY posterior_decoding(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *sb, PROBABILITY *sr, int obs_index, int state) {
  return sb[obs_index] * sr[obs_index] * fetch_forward_prob(model_def, sequence, f_table, obs_index, state) * fetch_backward_prob(model_def, sequence, b_table, obs_index, state);
//  return exp(log(sb[obs_index]) + sr[obs_index] + log(fetch_forward_prob(model_def, sequence, f_table, obs_index, state)) + log(fetch_backward_prob(model_def, sequence, b_table, obs_index, state)));
}


void calc_sr(PROBABILITY *sf, PROBABILITY *sb, int len, PROBABILITY *sr) {
  int i;

  sr[len - 1] = 1;
//  sr[len - 1] = log(1);
  for (i = len - 1; i > 0; i--) {
    sr[i - 1] = sr[i] * sb[i] / sf[i];
//    sr[i - 1] = sr[i] + log(sb[i]) - log(sf[i]);
  }
}


void print_forward_table(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, PROBABILITY *s, int n) {
  int i, j;
  PROBABILITY constant = 0;

  for (i = 0; i < sequence->len; i++) {
    constant += log(s[i]);
    if (i > n && i < sequence->len - n - 1) {
      fprintf(model_def->output, "%d\t", i);
      for (j = 0; j < model_def->n_states; j++) {
//        PROBABILITY p = log(fetch_forward_prob(model_def, sequence, table, i, j)) + constant;
//        p = exp(p);
        PROBABILITY p = fetch_forward_prob(model_def, sequence, table, i, j);
        if (j < model_def->n_states - 1) {
          fprintf(model_def->output, "%.14f\t", p);
        } else {
          fprintf(model_def->output, "%.14f\n", p);
        }
      }
    }
  }
}


void print_backward_table(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, PROBABILITY *s, int n) {
  int i, j;
  PROBABILITY constant = 0;

  for (i = sequence->len - 1; i >= 0; i--) {
    constant += log(s[i]);
    if (i > n && i < sequence->len - n - 1) {
      fprintf(model_def->output, "%d\t", (int)(sequence->len - 1 - i));
      for (j = 0; j < model_def->n_states; j++) {
//        PROBABILITY p = log(fetch_backward_prob(model_def, sequence, table, i, j)) + constant;
//        p = exp(p);
        PROBABILITY p = fetch_backward_prob(model_def, sequence, table, i, j);
        if (j < model_def->n_states - 1) {
          fprintf(model_def->output, "%.14f\t", p);
        } else {
          fprintf(model_def->output, "%.14f\n", p);
        }
      }
    }
  }
}


void print_posterior(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *s, int states_to_include, int n) {
  int i, j;

  for (i = 0; i < sequence->len; i++) {
    if (i < n || i > sequence->len - n - 1) {
      fprintf(model_def->output, "%d\t", i);
      for (j = 0; j < states_to_include; j++) {
        PROBABILITY p = log(s[i]) + log(fetch_forward_prob(model_def, sequence, f_table, i, j)) + log(fetch_backward_prob(model_def, sequence, b_table, i, j));
        if (j < states_to_include - 1) {
          fprintf(model_def->output, "%.14f\t", p);
        } else {
          fprintf(model_def->output, "%.14f\n", p);
        }
      }
    }
  }
}


void print_most_probable_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *s) {
  int i, j, max_index;
  PROBABILITY max_value;

  for (i = 0; i < sequence->len; i++) {
    fprintf(model_def->output, "%d\t", i);
    max_value = log(s[i]) + log(fetch_forward_prob(model_def, sequence, f_table, i, 0)) + log(fetch_backward_prob(model_def, sequence, b_table, i, 0));
    max_index = 0;
    for (j = 1; j < model_def->n_states; j++) {
      PROBABILITY p = log(s[i]) + log(fetch_forward_prob(model_def, sequence, f_table, i, j)) + log(fetch_backward_prob(model_def, sequence, b_table, i, j));
      if (p > max_value) {
        max_value = p;
        max_index = j;
      }
    }
    fprintf(model_def->output, "%d\n", max_index);
  }
}


PROBABILITY A_kl(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, int n_seqs, int k, int l) {
  int i, j;
  PROBABILITY sum = 0;
  PROBABILITY transition_prob = fetch_transition_prob(model_def, k, l);

  for (j = 0; j < n_seqs; j++) {
    for (i = 0; i < sequence[j]->len - 1; i++) {
      PROBABILITY emission_prob = fetch_emission_prob(model_def, l, sequence[j]->seq[i + 1]);
      PROBABILITY f_i = fetch_forward_prob(model_def, sequence[j], f_table[j], i, k);
      PROBABILITY b_i = fetch_backward_prob(model_def, sequence[j], b_table[j], i + 1, l);
      sum += f_i * transition_prob * emission_prob * b_i;
    }
  }

  return sum;
}


PROBABILITY E_kb(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, PROBABILITY **s, int n_seqs, int k, char b) {
  int i, j;
  PROBABILITY sum = 0;

  for (j = 0; j < n_seqs; j++) {
    for (i = 0; i < sequence[j]->len - 1; i++) {
      if (sequence[j]->seq[i] == b) {
        PROBABILITY f_i = fetch_forward_prob(model_def, sequence[j], f_table[j], i, k);
        PROBABILITY b_i = fetch_backward_prob(model_def, sequence[j], b_table[j], i, k);
        sum += s[j][i] * f_i * b_i;
      }
    }
  }

  return sum;
}


void print_transition_matrix_as_dot(model_def_struct *model_def) {
  int i, j;

  fprintf(model_def->output, "digraph G {\n  node [shape = circle];\n");
  for (i = 0; i < model_def->n_states; i++) {
    fprintf(model_def->output, "  node_%d [label = \"%d\"];\n", i, i);
  }

  fprintf(model_def->output, "\n");

  for (i = 0; i < model_def->n_states; i++) {
    for (j = 0; j < model_def->n_states; j++) {
      PROBABILITY p = fetch_transition_prob(model_def, i, j);
      if (p > 0.0) fprintf(model_def->output, "  node_%d -> node_%d [label = \"%.5f\"];\n", i, j, p);
    }
  }

  fprintf(model_def->output, "}\n");
}


void print_transition_matrix(model_def_struct *model_def) {
  int i, j;

  fprintf(model_def->output, "\n\t");
  for (i = 0; i < model_def->n_states - 1; i++) {
    fprintf(model_def->output, "%d\t", i);
  }
  fprintf(model_def->output, "%d\n", i);

  for (i = 0; i < model_def->n_states; i++) {
    fprintf(model_def->output, "%d\t", i);
    for (j = 0; j < model_def->n_states; j++) {
      fprintf(model_def->output, "%.5f", fetch_transition_prob(model_def, i, j));
      if (j < model_def->n_states - 1) {
        fprintf(model_def->output, "\t");
      } else {
        fprintf(model_def->output, "\n");
      }
    }
  }
}


void print_emission_matrix(model_def_struct *model_def) {
  int i, j;

  fprintf(model_def->output, "\n\t");
  for (i = 0; i < model_def->alphabet_length - 1; i++) {
    fprintf(model_def->output, "%d\t", i);
  }
  fprintf(model_def->output, "%d\n", i);

  for (i = 0; i < model_def->n_states; i++) {
    fprintf(model_def->output, "%d\t", i);
    for (j = 0; j < model_def->alphabet_length; j++) {
      fprintf(model_def->output, "%.5f", fetch_emission_prob(model_def, i, j));
      if (j < model_def->alphabet_length - 1) {
        fprintf(model_def->output, "\t");
      } else {
        fprintf(model_def->output, "\n");
      }
    }
  }
}


PROBABILITY log_likelihood(PROBABILITY **s, sequence_struct **sequence, int n_seqs) {
  int i, j;
  PROBABILITY sum = 0;

  for (i = 0; i < n_seqs; i++) {
    for (j = 0; j < sequence[i]->len; j++) {
      sum += log(s[i][j]);
    }
  }

  return sum;
}


void *forward_thread_wrapper(void *arg) {
  thread_wrapper_struct *tw = (thread_wrapper_struct *)arg;
  forward(tw->model_def, tw->sequence, tw->table, tw->s);
}


void *backward_thread_wrapper(void *arg) {
  thread_wrapper_struct *tw = (thread_wrapper_struct *)arg;
  backward(tw->model_def, tw->sequence, tw->s, tw->table);
}


void fb_on_all_seqs(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, PROBABILITY **sf, PROBABILITY **sb, int n_seqs) {
  int i, n_cpu;
  pthread_t *fpth, *bpth;
  thread_wrapper_struct *ftw, *btw;

  n_cpu = find_num_cpus();

  if (n_cpu < 2) {
    for (i = 0; i < n_seqs; i++) {
      forward(model_def, sequence[i], f_table[i], sf[i]);
      backward(model_def, sequence[i], sb[i], b_table[i]);
    }
  } else {
    fpth = ALLOC(sizeof(pthread_t) * n_seqs);
    bpth = ALLOC(sizeof(pthread_t) * n_seqs);
    ftw = ALLOC(sizeof(thread_wrapper_struct) * n_seqs);
    btw = ALLOC(sizeof(thread_wrapper_struct) * n_seqs);

    for (i = 0; i < n_seqs; i++) {
      ftw[i].model_def = model_def;
      ftw[i].sequence = sequence[i];
      ftw[i].table = f_table[i];
      ftw[i].s = sf[i];

      btw[i].model_def = model_def;
      btw[i].sequence = sequence[i];
      btw[i].table = b_table[i];
      btw[i].s = sb[i];

      pthread_create(bpth + i, NULL, backward_thread_wrapper, btw + i);
      pthread_create(fpth + i, NULL, forward_thread_wrapper, ftw + i);
    }

    for (i = 0; i < n_seqs; i++) {
      pthread_join(bpth[i], NULL);
      pthread_join(fpth[i], NULL);
    }

    free(fpth);
    free(bpth);
    free(ftw);
    free(btw);
  }
}


model_def_struct *initialize_model(char *filename, fixed_states_struct *fixed_states, int n_fixed_states) {
  model_def_struct *model_def;
  struct config_t cfg;
  config_setting_t *list;
  int i, j;

  config_init(&cfg);
  if (!config_read_file(&cfg, filename)) {
    fprintf(stderr, "Failed to load config file '%s'.\n", filename);
    exit(1);
  }

  model_def = ALLOC(sizeof(model_def_struct));
  model_def->n_fixed_states = 0;  // workaround of set_transition_prob looking for n_fixed_states to be set
  
  model_def->n_states = config_lookup_int(&cfg, "model.n_states");
//  fprintf(stderr, "n_states: %d\n", model_def->n_states);
  model_def->silent_states_begin = config_lookup_int(&cfg, "model.silent_states_begin");
  model_def->alphabet_length = config_lookup_int(&cfg, "model.alphabet_length");
  model_def->alphabet = ALLOC(sizeof(char) * model_def->alphabet_length);
  memcpy(model_def->alphabet, config_lookup_string(&cfg, "model.alphabet"), model_def->alphabet_length);

  model_def->initial_probs = ALLOC(sizeof(PROBABILITY) * model_def->silent_states_begin);
  list = config_lookup(&cfg, "model.initial_probs");
  for (i = 0; i < config_setting_length(list); i++) {
    config_setting_t * element = config_setting_get_elem(list, i);
    int state = config_setting_get_int_elem(element, 0);
    model_def->initial_probs[state] = config_setting_get_float_elem(element, 1);
  }

  model_def->transition_matrix = ALLOC(sizeof(PROBABILITY) * model_def->n_states * model_def->n_states);
  memset(model_def->transition_matrix, 0, sizeof(PROBABILITY) * model_def->n_states * model_def->n_states);
  list = config_lookup(&cfg, "model.transition_matrix");
  for (i = 0; i < config_setting_length(list); i++) {
    config_setting_t * element = config_setting_get_elem(list, i);
    int from = config_setting_get_int_elem(element, 0);
    int to = config_setting_get_int_elem(element, 1);
    set_transition_prob(model_def, from, to, config_setting_get_float_elem(element, 2));
  }
//  print_transition_matrix(model_def);

  if (n_fixed_states > 0) {
    model_def->fixed_states = fixed_states;
    model_def->n_fixed_states = n_fixed_states;
  } else {
    model_def->fixed_states = NULL;
    model_def->n_fixed_states = 0;
  }

  model_def->emission_matrix = ALLOC(sizeof(PROBABILITY) * model_def->alphabet_length * model_def->n_states);
  // the non-silent states default to 0
  memset(model_def->emission_matrix, 0, sizeof(PROBABILITY) * model_def->alphabet_length * model_def->silent_states_begin);
  list = config_lookup(&cfg, "model.emission_matrix");
  for (i = 0; i < config_setting_length(list); i++) {
    config_setting_t * element = config_setting_get_elem(list, i);
    int state = config_setting_get_int_elem(element, 0);
    int chr = config_setting_get_int_elem(element, 1);
    set_emission_prob(model_def, state, chr, config_setting_get_float_elem(element, 2));
  }

  // the silent states are all 1, but you can't memset this as a group..
  for (i = model_def->silent_states_begin; i < model_def->n_states; i++) {
    for (j = 0; j < model_def->alphabet_length; j++) {
      set_emission_prob(model_def, i, j, 1);
    }
  }
//  print_emission_matrix(model_def);


  model_def->parents = ALLOC(sizeof(int*) * model_def->n_states);
  model_def->n_parents = ALLOC(sizeof(int) * model_def->n_states);
  model_def->first_silent_parent = ALLOC(sizeof(int) * model_def->n_states);
  model_def->children = ALLOC(sizeof(int*) * model_def->n_states);
  model_def->n_children = ALLOC(sizeof(int) * model_def->n_states);
  model_def->first_silent_child = ALLOC(sizeof(int) * model_def->n_states);
  find_parents_and_children(model_def);

  config_destroy(&cfg);
  return model_def;
}


void find_parents_and_children(model_def_struct *model_def) {
  int *parents, *children;
  int n_parents, n_children, first_silent_parent, first_silent_child;
  int i, j;

  parents = ALLOC(sizeof(int) * model_def->n_states);
  children = ALLOC(sizeof(int) * model_def->n_states);

  for (i = 0; i < model_def->n_states; i++) {
    memset(parents, 0, sizeof(int) * model_def->n_states);
    memset(children, 0, sizeof(int) * model_def->n_states);
    n_parents = n_children = 0;

    for (j = 0; j < model_def->silent_states_begin; j++) {
      if (fetch_transition_prob(model_def, j, i) != 0.0) {
        parents[n_parents] = j;
        n_parents++;
      }
      if (fetch_transition_prob(model_def, i, j) != 0.0) {
        children[n_children] = j;
        n_children++;
      }
    }

    first_silent_parent = n_parents;
    first_silent_child = n_children;

    for (j = model_def->silent_states_begin; j < model_def->n_states; j++) {
      if (fetch_transition_prob(model_def, j, i) != 0.0) {
        parents[n_parents] = j;
        n_parents++;
      }
      if (fetch_transition_prob(model_def, i, j) != 0.0) {
        children[n_children] = j;
        n_children++;
      }
    }

    model_def->parents[i] = ALLOC(sizeof(int) * n_parents);
    memcpy(model_def->parents[i], parents, sizeof(int) * n_parents);
    model_def->n_parents[i] = n_parents;
    model_def->first_silent_parent[i] = first_silent_parent;

    model_def->children[i] = ALLOC(sizeof(int) * n_children);
    memcpy(model_def->children[i], children, sizeof(int) * n_children);
    model_def->n_children[i] = n_children;
    model_def->first_silent_child[i] = first_silent_child;
  }

  free(parents);
  free(children);
}


int read_sequence(char *filename, sequence_struct ***sequence_ptr) {
  sequence_struct **sequence;
  FILE *f, *f_index;
  int result, i, j, n_seqs = 0, begin_read, end_read;
  char **names, str[256];
  int max_filename_length = 256;

  f_index = fopen(filename, "r");
  // this is ugly, but I don't know how else to count the darn things
  while (fscanf(f_index, "%s %d %d\n", str, &begin_read, &end_read) > 0) {
//    fprintf(stderr, "file %d: \"%s\"\n", n_seqs, str);
    n_seqs++;
  }
  rewind(f_index);

  sequence = ALLOC(sizeof(sequence_struct *) * n_seqs);

  i = 0;
  while (fscanf(f_index, "%s %d %d\n", str, &begin_read, &end_read) > 0) {
    sequence[i] = ALLOC(sizeof(sequence_struct));


    if (end_read < 0) {
      struct stat sb;
      result = stat(str, &sb);
      if (result != 0) {
        fprintf(stderr, "Error reading %s.  Exiting.\n", str);
        exit;
      }

      end_read = sb.st_size;
    }

    sequence[i]->len = end_read - begin_read + 1;
    sequence[i]->seq = ALLOC(sizeof(char) * sequence[i]->len);

    f = fopen(str, "r");
    fseek(f, begin_read - 1, SEEK_SET);
    result = fread(sequence[i]->seq, sequence[i]->len, 1, f);
//    fprintf(stderr, "%d objects of size %d read from %s.\n", result, sequence[i]->len, str);
    fclose(f);

    i++;
  }

  fclose(f_index);

  *sequence_ptr = sequence;
  return n_seqs;
}


void print_initial_probs(model_def_struct *model_def) {
  int i;

  fprintf(model_def->output, "\ninitial_probs:\n");
  for (i = 0; i < model_def->silent_states_begin; i++) {
    fprintf(model_def->output, "%d:\t%f\n", i, model_def->initial_probs[i]);
  }
}


BOOL verify_model(model_def_struct *model_def) {
  int i, j;
  BOOL error = FALSE;
  PROBABILITY t_sum, e_sum;
  PROBABILITY eps = 1e-5;

  for (i = 0; i < model_def->n_states; i++) {
    t_sum = 0.0;
    e_sum = 0.0;

    for (j = 0; j < model_def->n_states; j++) {
      t_sum += fetch_transition_prob(model_def, i, j);
    }
    for (j = 0; j < model_def->alphabet_length; j++) {
      e_sum += fetch_emission_prob(model_def, i, j);
    }

    if (fabs(t_sum - 1.0) > eps) {
      printf ("Transitions from state %d sum to %f instead of 1.\n", i, t_sum);
      error = TRUE;
    }
    if (fabs(e_sum - 1.0) > eps) {
      printf ("Emissions from state %d sum to %f instead of 1.\n", i, e_sum);
      error = TRUE;
    }
  }

  return error;
}


int find_num_cpus() {
  int nprocs = -1;
#ifdef _SC_NPROCESSORS_ONLN
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);
#elif defined CTL_HW
  int mib[2];
  size_t len;

  mib[0] = CTL_HW;
  mib[1] = HW_AVAILCPU;
  len = sizeof(nprocs);
  sysctl(mib, 2, &nprocs, &len, NULL, 0);
#endif

  return nprocs;
}


static int cmp_sortable_pairs(const void *a, const void *b) {
  if (((sortable_pair *)a)->value > ((sortable_pair *)b)->value) return 1;
  else if (((sortable_pair *)a)->value < ((sortable_pair *)b)->value) return -1;
  else return 0;
}


void update_silent_states_viterbi(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int row_index, sortable_pair *pairs, int *ptr) {
  PROBABILITY *row = table + (unsigned long)model_def->n_states * (unsigned long)row_index;
  int *argmax = ptr + (unsigned long)model_def->n_states * (unsigned long)row_index;
  int i, j, n_pairs;

  for (i = model_def->silent_states_begin; i < model_def->n_states; i++) {
    n_pairs = model_def->first_silent_parent[i];
    for (j = 0; j < n_pairs; j++) {
      int state = model_def->parents[i][j];
      PROBABILITY transition_prob = fetch_transition_prob(model_def, state, i);
      pairs[j].index = state;
      pairs[j].value = row[state] * transition_prob;
    }

    qsort(pairs, n_pairs, sizeof(sortable_pair), cmp_sortable_pairs);
    row[i] = pairs[n_pairs-1].value;
    argmax[i] = pairs[n_pairs-1].index;
  }

  for (i = model_def->silent_states_begin + 1; i < model_def->n_states; i++) {
    n_pairs = model_def->n_parents[i] - model_def->first_silent_parent[i];
    for (j = model_def->first_silent_parent[i]; j < model_def->n_parents[i]; j++) {
      int state = model_def->parents[i][j];
      PROBABILITY transition_prob = fetch_transition_prob(model_def, state, i);
      pairs[j].index = state;
      pairs[j].value = row[state] * transition_prob;
    }

    if (n_pairs > 0) {
      qsort(pairs, n_pairs, sizeof(sortable_pair), cmp_sortable_pairs);
      if (pairs[n_pairs-1].value > row[i]) {
        row[i] = pairs[n_pairs-1].value;
        argmax[i] = pairs[n_pairs-1].index;
      }
    }
  }
}


void update_normal_states_viterbi(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int row_index, sortable_pair *pairs, int *ptr) {
  PROBABILITY *row = table + (unsigned long)model_def->n_states * (unsigned long)row_index;
  PROBABILITY *prev_row = table + (unsigned long)model_def->n_states * (unsigned long)(row_index - 1);
  int *argmax = ptr + (unsigned long)model_def->n_states * (unsigned long)row_index;
  int i, j, n_pairs;

  for (i = 0; i < model_def->silent_states_begin; i++) {
    n_pairs = model_def->n_parents[i];
    for (j = 0; j < n_pairs; j++) {
      int state = model_def->parents[i][j];
      PROBABILITY transition_prob = fetch_transition_prob(model_def, state, i);
      pairs[j].index = state;
      pairs[j].value = prev_row[state] * transition_prob;
    }

    qsort(pairs, n_pairs, sizeof(sortable_pair), cmp_sortable_pairs);
    row[i] = fetch_emission_prob(model_def, i, sequence->seq[row_index]) * pairs[n_pairs-1].value;
    argmax[i] = pairs[n_pairs-1].index;
  }
}

void viterbi(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *table, int *path, PROBABILITY *s) {
  int i, j;
  sortable_pair *pairs;
  int *ptr;

  pairs = ALLOC(sizeof(sortable_pair) * model_def->n_states);
  ptr = ALLOC(sizeof(int) * model_def->n_states * sequence->len);

  // initialize first row
  // normal states need to be handled specially
  for (i = 0; i < model_def->silent_states_begin; i++) {
    table[i] = model_def->initial_probs[i] * fetch_emission_prob(model_def, i, sequence->seq[0]);
    ptr[i] = model_def->n_states;  // special value to indicate no previous state
  }

  // silent states can use the normal machinery
  update_silent_states_viterbi(model_def, sequence, table, 0, pairs, ptr);

  // the first weight is just the sum of the first column (row, in this implementation).  calculate, and normalize.
  s[0] = normalize_row(table, model_def->silent_states_begin, model_def->n_states);

  // first row initialization complete


  // fill out the rest of the table now
  for (i = 1; i < sequence->len; i++) {
    #ifdef VERBOSE
    if (i % 500 == 0) fprintf(stderr, "viterbi row %d\n", i);
    #endif
    update_normal_states_viterbi(model_def, sequence, table, i, pairs, ptr);
    update_silent_states_viterbi(model_def, sequence, table, i, pairs, ptr);
    s[i] = normalize_row(table + (unsigned long)model_def->n_states * (unsigned long)i, model_def->silent_states_begin, model_def->n_states);
  }

  PROBABILITY *row = table + (unsigned long)model_def->n_states * (unsigned long)(sequence->len-1);
  int state;
  for (i = 0; i < model_def->n_states; i++) {
    pairs[i].index = i;
    pairs[i].value = row[i];
  }
  qsort(pairs, model_def->n_states, sizeof(sortable_pair), cmp_sortable_pairs);
  state = pairs[model_def->n_states-1].index;

  for (i = sequence->len - 1; i >= 0; i--) {
    int *argmax = ptr + (unsigned long)model_def->n_states * (unsigned long)(i);
    while (state >= model_def->silent_states_begin) { // backtrack through silent states
      state = argmax[state];
    }
    path[i] = state;
    state = argmax[state];
  }

  free(pairs);
  free(ptr);
}



void load_seq_pos_conc_scaler(char* file_name, float** seq_pos_conc_scaler,
                              int sequence_length, int n_motifs, BOOL nuc_present)
{
    // the file that specifies the position specific conc scaler should be a csv
    // file, delimited by \t. The first line of the file should be header.
    // The first column should be nucleosome, and the rest should be TFs
    FILE *in;
    //char tmp[1024];
    char tmp[4096];

    int line_count = 0;
    int field_parsed = 0;

    in = fopen(file_name, "r");

    if(in == NULL)
    {
        fprintf(stderr, "Can't open %s\n", file_name);
        exit(1);
    }

    // ignore the header line
    char* ptr = fgets(tmp, sizeof(tmp),in);
    line_count++;

    // how many filed do I expect on each line
    int total_filed = n_motifs;
    if(nuc_present)
    {
        total_filed = n_motifs + 1;
    }

    while(fgets(tmp, sizeof(tmp),in))
    {
        if(line_count > sequence_length)
        {
            fprintf(stderr, "hit line %d in file %s, but sequence length is %d\n",
                    line_count, file_name, sequence_length);
            exit(1);
        }

        field_parsed = parse_one_line(tmp, seq_pos_conc_scaler, line_count - 1,
                                      total_filed);

        // filed in each line should equal to the number of TFs plus nucleosome
        //if(field_parsed != total_filed) {
        //  fprintf(stderr, "line %d only have %d field, there should be %d\n",
        //          line_count, field_parsed, n_motifs + 1);
        //  exit(1);
        //}

        line_count++;
    }

    if((line_count - 1) != sequence_length)
    {
        fprintf(stderr, "file %s only has %d lines, but sequence length is %d\n",
                file_name, line_count, sequence_length);
        exit(1);
    }

    fclose(in);
}

int parse_one_line(char* line, float** matrix, int row, int total_fields)
{
    int field_parsed = 0;
    //char tmp[10];

    char *i, *j;

    j = line;
    i = line;
    while(*j != '\0')
    {
        if(*j == '\t' || *j == '\n')
        {
            char tmp[j-i+1];
            strncpy(tmp, i, (j - i));
            // without this '\0', there would be a bug
            tmp[j-i] = '\0';
            if(field_parsed > total_fields)
            {
                fprintf(stderr, "filed number exceeds total field, %d > %d at line %d\n",
                        field_parsed, total_fields, row + 1);
                exit(1);

            }
            matrix[field_parsed][row] = atof(tmp);
            //fprintf(stderr,
            //              "TEST: At row %ld, for filed_parsed %ld, the scaling factor is %s as string and  %.12f as value\n",
            //              row, field_parsed, tmp, matrix[field_parsed][row]);
            field_parsed++;
            i = j + 1;
        }
        j++;
    }

    return field_parsed;
}