#include "bc.h"
#include <time.h>
#include <libgen.h>
#include <getopt.h>

extern char *optarg;
extern int optind;

void free_memory(model_def_struct *model_def, sequence_struct **sequence, PROBABILITY **f_table, PROBABILITY **b_table, PROBABILITY **sf, PROBABILITY **sb, PROBABILITY **sr, int n_seqs, int *motif_starts, int *motif_lens, PROBABILITY *motif_conc, int n_motifs, char **motif_names) {
  int i;

  free(model_def->initial_probs);
  free(model_def->transition_matrix);
  if (model_def->n_fixed_states > 0) {
    free(model_def->fixed_states);
  }
  free(model_def->emission_matrix);
  free(model_def->alphabet);
  free(model_def->n_parents);
  free(model_def->n_children);
  free(model_def->first_silent_parent);
  free(model_def->first_silent_child);
  for (i = 0; i < model_def->n_states; i++) {
    free(model_def->parents[i]);
    free(model_def->children[i]);
  }
  free(model_def->parents);
  free(model_def->children);
  free(model_def);

  for (i = 0; i < n_seqs; i++) {
    free(sequence[i]->seq);
    free(sequence[i]);
    free(f_table[i]);
    free(b_table[i]);
    free(sf[i]);
    free(sb[i]);
    free(sr[i]);
  }

  if (motif_names != NULL) {
    for (i = 0; i < n_motifs; i++) {
      free(motif_names[i]);
    }
  }
  free(motif_names);

  free(sequence);
  free(f_table);
  free(b_table);
  free(motif_starts);
  free(motif_lens);
  free(sf);
  free(sb);
  free(sr);
  free(motif_conc);
}


void find_motif_state_numbers(model_def_struct *model_def, int **starts, int **lens) {
  int *motif_starts, *motif_lens;
  int i, j;

  motif_starts = ALLOC(sizeof(int) * (model_def->n_states - model_def->silent_states_begin));
  motif_lens = ALLOC(sizeof(int) * (model_def->n_states - model_def->silent_states_begin));

  for (i = model_def->silent_states_begin + 1; i < model_def->n_states; i++) {
    // i is each silent state for a motif
    int first = -1, second = -1;
    for (j = 1; j < model_def->silent_states_begin; j++) {
      if (fetch_transition_prob(model_def, i, j) > 0) {
        if (first == -1) {
          first = j;
        } else {
          second = j;
          break;
        }
      }
    }

    motif_starts[i - model_def->silent_states_begin - 1] = first > second ? second : first;
    motif_lens[i - model_def->silent_states_begin - 1] = abs(first - second);
//    fprintf(stderr, "Motif %d start: %d\tlen: %d\n", i, motif_starts[i - model_def->silent_states_begin - 1], motif_lens[i - model_def->silent_states_begin - 1]);
  }

  *starts = motif_starts;
  *lens = motif_lens;
}


BOOL find_nucleosome_states(model_def_struct *model_def, int *motif_starts, int *motif_lens, int *nuc_start, int *nuc_len) {
  int distributor_index = model_def->silent_states_begin;
  int n_motifs = model_def->n_states - model_def->silent_states_begin - 1;
  int end_of_last_motif = 1;
  if (n_motifs > 0) end_of_last_motif = motif_starts[n_motifs-1] + 2*motif_lens[n_motifs-1];

  if ((*nuc_len = distributor_index - end_of_last_motif) > 0) {
    *nuc_start = distributor_index - *nuc_len;
    return TRUE;
  } else {
    return FALSE;
  }
}


int find_num_nucleosome_padding_states(model_def_struct *model_def, int nuc_start) {
// the background states have transitions of exactly 1 leading out of them and into the next state.
// following this along will show how many there are, plus 4 for the branched background state.

  int i = nuc_start;
  while (fetch_transition_prob(model_def, i, i + 1) == 1.0) i++;

  return (i - nuc_start + 5);
}


void update_a0k_probabilities(model_def_struct *model_def) {
// I know the silent states before each motif start at model_def->silent_states_begin + 1
// I also know that the transitions from each one of these goes to exactly two states, which are the beginning of the forward and reverse motifs respectively.  So abs(index_of_one - index_of_other) must be the length of the motifs.
// I also know the transition probabilities from the distributor state into each of these.
// Therefore I have the information I need to implement the updating of a_{0k} as described in my document.

  int *motif_starts, *motif_lens;
  int n_motifs = model_def->n_states - model_def->silent_states_begin - 1;
  int distributor_index = model_def->silent_states_begin;
  int i, j;
  PROBABILITY sum = 0;

  find_motif_state_numbers(model_def, &motif_starts, &motif_lens);
  // start with background component, which is only 1 long...
  sum = fetch_transition_prob(model_def, distributor_index, 0);
  model_def->initial_probs[0] = sum;

  for (i = 0; i < n_motifs; i++) {
    PROBABILITY p = fetch_transition_prob(model_def, distributor_index, distributor_index + i + 1);
    for (j = motif_starts[i]; j < motif_starts[i] + 2 * motif_lens[i]; j++) {
      model_def->initial_probs[j] = p;
      sum += p;
    }
  }

  // determine if nucleosome is present
  int nuc_len = 0;
  int nuc_start = 0;
  if (find_nucleosome_states(model_def, motif_starts, motif_lens, &nuc_start, &nuc_len)) {
    PROBABILITY p = fetch_transition_prob(model_def, distributor_index, nuc_start);
    int n_padding_states = find_num_nucleosome_padding_states(model_def, nuc_start);

    // left (normal) padding states
    for (i = nuc_start; i < nuc_start + n_padding_states - 4; i++) {
      model_def->initial_probs[i] = p;
      sum += p;
    }

    // branched padding state
    for (i = nuc_start + n_padding_states - 4; i < nuc_start + n_padding_states; i++) {
      model_def->initial_probs[i] = p / 4.0;
      sum += p / 4.0;
    }

    // tons of nucleosome states, 16 per sequence position
    for (i = nuc_start + n_padding_states; i < nuc_start + nuc_len - n_padding_states + 3; i++) {
      model_def->initial_probs[i] = p / 16.0;
      sum += p / 16.0;
    }

    // right branching states, all of which are normal
    for (i = nuc_start + nuc_len - n_padding_states + 3; i < nuc_start + nuc_len; i++) {
      model_def->initial_probs[i] = p;
      sum += p;
    }
  }

  for (i = 0; i < model_def->silent_states_begin; i++) {
    model_def->initial_probs[i] /= sum;
  }

  free(motif_starts);
  free(motif_lens);
}


void apply_temperature(model_def_struct *model_def, int *motif_starts, int *motif_lens, int nuc_start, int nuc_len, PROBABILITY T) {
// T is the inverse temperature parameter, as describe in Segal's ImplementationNotes.pdf
  int i, j, k;

  // scale background emissions
  for (i = 0; i < model_def->alphabet_length; i++) {
    PROBABILITY p = fetch_emission_prob(model_def, 0, i);
    set_emission_prob(model_def, 0, i, pow(p, T));
  }

  // scale motif emissions
  int n_motifs = model_def->n_states - model_def->silent_states_begin - 1;
  for (i = 0; i < n_motifs; i++) {
    for (j = motif_starts[i]; j < motif_starts[i] + 2 * motif_lens[i]; j++) {
      for (k = 0; k < model_def->alphabet_length; k++) {
        PROBABILITY p = fetch_emission_prob(model_def, j, k);
        set_emission_prob(model_def, j, k, pow(p, T));
      }
    }
  }


  // scale nucleosome emissions/transitions
  if (nuc_start > 0) {
    int n_padding_states = find_num_nucleosome_padding_states(model_def, nuc_start);

    // number of bases/positions in the actual nucleosome; ideally 147 but practically less due to data limitations; padding on the right is 3 shorter than the left, since there's no branched bg state
    int n_nuc_pos = (nuc_len - (2 * n_padding_states - 3)) / 16;

    // normal background states in the beginning of the padding
    for (i = nuc_start; i < nuc_start + n_padding_states - 4; i++) {
      for (j = 0; j < model_def->alphabet_length; j++) {
        PROBABILITY p = fetch_emission_prob(model_def, i, j);
        set_emission_prob(model_def, i, j, pow(p, T));
      }
    }

    // normal background states at the end of the padding
    for (i = nuc_start + nuc_len - (n_padding_states - 3); i < nuc_start + nuc_len; i++) {
      for (j = 0; j < model_def->alphabet_length; j++) {
        PROBABILITY p = fetch_emission_prob(model_def, i, j);
        set_emission_prob(model_def, i, j, pow(p, T));
      }
    }

    // transitions into branched background state's 4 branches
    for (i = nuc_start + n_padding_states - 4; i < nuc_start + n_padding_states; i++) {
      PROBABILITY p = fetch_transition_prob(model_def, nuc_start + n_padding_states - 5, i);
      set_transition_prob(model_def, nuc_start + n_padding_states - 5, i, pow(p, T));
    }

    // handle the transitions from the branched background state into the first nucleosome states
    for (i = nuc_start + n_padding_states - 4; i < nuc_start + n_padding_states; i++) {
      for (j = nuc_start + n_padding_states; j < nuc_start + n_padding_states + 16; j++) {
        PROBABILITY p = fetch_transition_prob(model_def, i, j);
        set_transition_prob(model_def, i, j, pow(p, T));
      }
    }

    // handle the transitions between nucleosome states
    for (i = 0; i < n_nuc_pos - 1; i++) {
      for (j = nuc_start + n_padding_states + 16 * i; j < nuc_start + n_padding_states + 16 * (i + 1); j++) {
        for (k = nuc_start + n_padding_states + 16 * (i + 1); k < nuc_start + n_padding_states + 16 * (i + 2); k++) {
          PROBABILITY p = fetch_transition_prob(model_def, j, k);
          set_transition_prob(model_def, j, k, pow(p, T));
        }
      }
    }

  }
}


PROBABILITY posterior_over_range(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *sb, PROBABILITY *sr, int pos, int from, int to) {
  int i;
  PROBABILITY sum = 0;

  for (i = from; i <= to; i++) {
    sum += posterior_decoding(model_def, sequence, f_table, b_table, sb, sr, pos, i);
  }

  return sum;
}


void posterior_output_summed_states(model_def_struct *model_def, sequence_struct *sequence, PROBABILITY *f_table, PROBABILITY *b_table, PROBABILITY *sb, PROBABILITY *sr, int *motif_starts, int *motif_lens, char **motif_names, BOOL output_start_probs_only) {
  int i, j, k;
  int n_motifs = model_def->n_states - model_def->silent_states_begin - 1;
  int nuc_start, nuc_len, n_padding_states;
  BOOL nuc_present = FALSE;

  if (nuc_present = find_nucleosome_states(model_def, motif_starts, motif_lens, &nuc_start, &nuc_len)) {
//    n_padding_states = find_num_nucleosome_padding_states(model_def, nuc_start);
    n_padding_states = 5;
  }

  // print header
  fprintf(model_def->output, "background");
  if (!motif_names) {
    for (j = 0; j < n_motifs; j++) {
      fprintf(model_def->output, "\tmotif_%d", j);
    }
  } else {
    for (j = 0; j < n_motifs; j++) {
      fprintf(model_def->output, "\t%s", motif_names[j]);
    }
  }
  if (nuc_present) fprintf(model_def->output, "\tnuc_padding\tnucleosome");
  fprintf(model_def->output, "\n");

  for (i = 0; i < sequence->len; i++) {
    fprintf(model_def->output, "%.20f\t", posterior_decoding(model_def, sequence, f_table, b_table, sb, sr, i, 0));

    for (j = 0; j < n_motifs; j++) {
      PROBABILITY sum = 0;
      if (!output_start_probs_only) {
        sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, motif_starts[j], motif_starts[j] + motif_lens[j] - 1);
        sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, motif_starts[j] + motif_lens[j], motif_starts[j] + 2 * motif_lens[j] - 1);
      } else {
        sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, motif_starts[j], motif_starts[j]);
        sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, motif_starts[j] + motif_lens[j], motif_starts[j] + motif_lens[j]);
      }
      if (j < n_motifs - 1) {
        fprintf(model_def->output, "%.20f\t", sum);
      } else {
        fprintf(model_def->output, "%.20f", sum);
      }
    }

    if (nuc_present) {
      PROBABILITY sum = 0;

      // calculate nuc_padding prob
      sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, nuc_start, nuc_start + n_padding_states - 1);
      sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, nuc_start + nuc_len - n_padding_states, nuc_start + nuc_len - 1);
      fprintf(model_def->output, "\t%.20f", sum);

      // calculate nucleosome prob
      sum = 0;
      if (!output_start_probs_only) {
        sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, nuc_start + n_padding_states, nuc_start + nuc_len - n_padding_states - 1);
      } else {
        sum += posterior_over_range(model_def, sequence, f_table, b_table, sb, sr, i, nuc_start + n_padding_states, nuc_start + n_padding_states);
      }
      fprintf(model_def->output, "\t%.20f", sum);
    }

    fprintf(model_def->output, "\n");
  }
}


void print_usage(char **argv) {
  fprintf(stderr, "usage: %s [options] model_file seq_file\n", basename(argv[0]));
  fprintf(stderr, "  -n  nucleosome_concentration (float)\n");
  fprintf(stderr, "  -m  motif_concentrations (comma delimited string of floats)\n");
  fprintf(stderr, "  -N  motif_labels (comma delimited string of strings, for output file column headers)\n");
  fprintf(stderr, "  -u  unbound_concentration (float)\n");
  fprintf(stderr, "  -t  inverse_temperature (float)\n");
  fprintf(stderr, "  -s  output only probabilities of starting each DBF per postion\n");
  fprintf(stderr, "\nexample: %s -n 1.0 -m 0.01,0.1,0.01 -u 1.0 -t 2.0 model.cfg seq_filenames.txt > output.txt\n", basename(argv[0]));
}


void parse_opts(int argc, char **argv, PROBABILITY *nuc_conc, PROBABILITY *unbound_conc, PROBABILITY *motif_conc, PROBABILITY *T, char **motif_names, char *fixed_states_str, BOOL *output_start_probs_only) {
  int opt, i;
  char *str, *token;

  while ((opt = getopt(argc, argv, "n:m:u:t:hN:s")) > 0) {
    switch (opt) {
      case 'n':
        *nuc_conc = atof(optarg);
        break;
      case 'u':
        *unbound_conc = atof(optarg);
        break;
      case 't':
        *T = atof(optarg);
        break;
      case 'm':
        for (i = 0, str = optarg; ; i++, str = NULL) {
          token = strtok(str, ",");
          if (token == NULL) break;
          motif_conc[i] = atof(token);
        }
        break;
      case 'N':
        for (i = 0, str = optarg; ; i++, str = NULL) {
          token = strtok(str, ",");
          if (token == NULL) break;
          motif_names[i] = ALLOC(strlen(token) * sizeof(char));
          strcpy(motif_names[i], token);
        }
        break;
      case 's':
        *output_start_probs_only = TRUE;
        break;

      case '?':
      case 'h':
      default:
        print_usage(argv);
        exit(0);
    }
  }
}


int main(int argc, char **argv) {
  model_def_struct *model_def;
  sequence_struct **sequence;
  int n_seqs = 0, i, j;
  PROBABILITY **f_table, **b_table, **sf, **sb, **sr;
  PROBABILITY T = 1.0;
  PROBABILITY nuc_conc = 1.0, unbound_conc = 1.0, *motif_conc;
  char **motif_names;
  BOOL nuc_present = FALSE;
  int n_fixed_states;
  int *fixed_states;


  if (argc < 3) {
    print_usage(argv);
    return 0;
  }

  motif_conc = ALLOC(sizeof(PROBABILITY) * 256); // I need to know the number of motifs before parsing, but I can't (easily), so for now I guess a max
  for (i = 0; i < 256; i++) motif_conc[i] = 0.01;
  motif_names = ALLOC(sizeof(char*) * 256); // I need to know the number of motifs before parsing, but I can't (easily), so for now I guess a max
  memset(motif_names, 0, sizeof(char*) * 256);

  // I'll leave the fixed state logic in here under the hood, but I've removed the interface to it to reduce confusion in these releases
  char fixed_states_str[256] = {'\0'};
  BOOL output_start_probs_only = FALSE;
  parse_opts(argc, argv, &nuc_conc, &unbound_conc, motif_conc, &T, motif_names, fixed_states_str, &output_start_probs_only);

  if (motif_names[0] == 0) {  // if -N wasn't on the command line, free this up so the output routine doesn't try to use it later
    free(motif_names);
    motif_names = NULL;
  }

  model_def = initialize_model(argv[optind], NULL, 0);
  n_seqs = read_sequence(argv[optind + 1], &sequence);

  f_table = ALLOC(sizeof(PROBABILITY *) * n_seqs);
  b_table = ALLOC(sizeof(PROBABILITY *) * n_seqs);
  sf = ALLOC(sizeof(PROBABILITY *) * n_seqs);
  sb = ALLOC(sizeof(PROBABILITY *) * n_seqs);
  sr = ALLOC(sizeof(PROBABILITY *) * n_seqs);
  for (i = 0; i < n_seqs; i++) {
    f_table[i] = ALLOC(sizeof(PROBABILITY) * model_def->n_states * sequence[i]->len);
    b_table[i] = ALLOC(sizeof(PROBABILITY) * model_def->n_states * sequence[i]->len);
    sf[i] = ALLOC(sizeof(PROBABILITY) * sequence[i]->len);
    sb[i] = ALLOC(sizeof(PROBABILITY) * sequence[i]->len);
    sr[i] = ALLOC(sizeof(PROBABILITY) * sequence[i]->len);
    memset(sf[i], 0, sizeof(PROBABILITY) * sequence[i]->len);
    memset(sb[i], 0, sizeof(PROBABILITY) * sequence[i]->len);
    memset(sr[i], 0, sizeof(PROBABILITY) * sequence[i]->len);
  }


//  if (!verify_model(model_def)) return -1;
  int nuc_len = 0;
  int nuc_start = 0;
  int n_motifs = model_def->n_states - model_def->silent_states_begin - 1;
  int states_len = 2 * n_motifs + 1;
  int *motif_starts, *motif_lens;
  find_motif_state_numbers(model_def, &motif_starts, &motif_lens);
  if (nuc_present = find_nucleosome_states(model_def, motif_starts, motif_lens, &nuc_start, &nuc_len)) states_len++;

  fprintf(stderr, "Inverse Temp.: %f\n", T);
  fprintf(stderr, "Unbound Conc.: %f\n", unbound_conc);
  if (nuc_present) fprintf(stderr, "Nucleosome Conc.: %f\n", nuc_conc);
  for (i = 0; i < n_motifs; i++)
    fprintf(stderr, "Motif %d Conc.: %f\n", i, motif_conc[i]);

  if (argc - optind > 2) {
    if (!(model_def->output = fopen(argv[optind + 2], "w"))) {
      fprintf(stderr, "Opening %s for writing failed.\n", argv[optind + 2]);
      exit(0);
    }
  } else {
    model_def->output = stdout;
  }

  set_transition_prob(model_def, model_def->silent_states_begin, 0, unbound_conc);
  set_transition_prob(model_def, model_def->silent_states_begin, nuc_start, nuc_conc);
  for (i = 0; i < n_motifs; i++)
    set_transition_prob(model_def, model_def->silent_states_begin, model_def->silent_states_begin + i + 1, motif_conc[i]);

  apply_temperature(model_def, motif_starts, motif_lens, nuc_start, nuc_len, T);
  update_a0k_probabilities(model_def);


  fprintf(stderr, "Running forward/backward and posterior decoding.\n");
  fb_on_all_seqs(model_def, sequence, f_table, b_table, sf, sb, n_seqs);

  calc_sr(sf[0], sb[0], sequence[0]->len, sr[0]);
  posterior_output_summed_states(model_def, sequence[0], f_table[0], b_table[0], sb[0], sr[0], motif_starts, motif_lens, motif_names, output_start_probs_only);

  fclose(model_def->output);
  free_memory(model_def, sequence, f_table, b_table, sf, sb, sr, n_seqs, motif_starts, motif_lens, motif_conc, n_motifs, motif_names);

  return 0;
}
