#include "bc.h"


void free_memory(model_def_struct *model_def) {
  free(model_def->initial_probs);
  free(model_def->transition_matrix);
  free(model_def->emission_matrix);
  free(model_def);
}


int main(int argc, char **argv) {
  model_def_struct *model_def;

  if (argc < 2) {
    printf("Usage: %s model_file\n", argv[0]);
    return 1;
  }

  model_def = initialize_model(argv[1], NULL, 0);
  model_def->output = stdout;
  print_transition_matrix_as_dot(model_def);

  free_memory(model_def);

  return 0;
}
