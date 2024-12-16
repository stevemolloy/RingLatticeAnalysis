#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "lib.h"
#include "elements.h"

#define M_PI (acos(-1.0))
#define ANGLE_EPSILON (0.1*M_PI/180)

#include "stb_ds.h"

#include "sdm_lib.h"

// TODO: Factor out the twiss propagation stuff into diagonalised 2x2 instead of the 3x3 stuff
// TODO: Rationalise memory management. In some places it is dynamic arrays, in other an arena is used. Time to straighten this out.
// TODO: Spin off a few more functions for the library.  The "main" function should be a list of calls to library functions.
// TODO: Rationalise organisation of the libraries.
// TODO: Handle off-energy dynamics.

void usage(const char* program) {
  fprintf(stderr, "USAGE: %s [-p <periodicity>] [-h <harmonic number>] [-E <kinetic energy>] lattice_file\n", program);
}

bool lattice_is_closed(double total_angle) {
  return !(fabs(total_angle - 2*M_PI) > ANGLE_EPSILON);
}

int main(int argc, char **argv) {
  CommandLineArgs args = get_clargs(argc, argv);
  if (args.error) {
    usage(args.programname);
    return 1;
  }

  const double gamma_0 = args.E_0 * 1e9 / ELECTRON_MASS;

  Element *line = {0};
  generate_lattice(args.file_path, &line);

  double line_length = calculate_line_length(line);
  double total_length = line_length * args.periodicity;

  double line_angle = calculate_line_angle(line);
  double total_angle = line_angle * args.periodicity;
  bool closed_system = lattice_is_closed(total_angle);

  printf("\nSummary of the lattice defined in %s\n\n", args.file_path);

  if (!closed_system) {
    printf("Total bending angle (%0.3f degrees) is not within %0.1e degrees of a circle.\n", 
           radians_to_degrees(total_angle), radians_to_degrees(ANGLE_EPSILON));
    printf("System does not close, so not calculating ring parameters.\n\n");
  }

  double line_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  double x_trace, y_trace, I[5] = {0};
  double I_1, I_2, I_3, I_4, I_5;
  double j_x, T_0;
  LinOptsParams lin_opt_params = {0};
  if (closed_system) {
    propagate_linear_optics(line, line_matrix, &lin_opt_params, I);
    apply_matrix_n_times(total_matrix, line_matrix, args.periodicity);

    x_trace = (total_matrix[0*BEAM_DOFS + 0] + total_matrix[1*BEAM_DOFS + 1]);
    y_trace = (total_matrix[2*BEAM_DOFS + 2] + total_matrix[3*BEAM_DOFS + 3]);

    if (args.save_twiss) {
      FILE *twiss_file = fopen(args.twiss_filename, "w");
      fprintf(twiss_file, "S / m, beta_x / m, beta_y / m, eta_x / m\n");
      for (size_t i=0; i<arrlenu(line); i++) {
        fprintf(twiss_file, "%0.6e, %0.6e, %0.6e, %0.6e\n",
		    		lin_opt_params.Ss[i],
		    		lin_opt_params.element_beta_xs[i],
		    		lin_opt_params.element_beta_ys[i],
		    		lin_opt_params.element_etas[i]);
      }
      fclose(twiss_file);
    }
    I_1 = I[0] * args.periodicity;
    I_2 = I[1] * args.periodicity;
    I_3 = I[2] * args.periodicity;
    I_4 = I[3] * args.periodicity;
    I_5 = I[4] * args.periodicity;

    j_x = 1.0f - I_4/I_2;
    T_0 = total_length / C;

    arrfree(lin_opt_params.element_etas);
    arrfree(lin_opt_params.element_etaps);
    arrfree(lin_opt_params.element_beta_xs);
    arrfree(lin_opt_params.element_beta_ys);
    arrfree(lin_opt_params.element_curlyH);
    arrfree(lin_opt_params.Ss);
  } else {
    get_line_matrix(line_matrix, line);
    apply_matrix_n_times(total_matrix, line_matrix, args.periodicity);
  }

  if (args.E_0 == 0) 
    printf("WARNING: Kinetic energy not provided, so not calculating all parameters\n\n");

  printf("Periodicity: %zu\n", args.periodicity);
  printf("Number of elements in the line: %td\n", arrlen(line));
  if (args.periodicity != 1) {
    printf("Total length of the lattice: %0.3f m (%0.3f m for the line)\n", total_length, line_length);
  } else {
    printf("Total length of the line: %f m\n", line_length);
  }
  if (args.periodicity != 1) {
    printf("Total bending angle of the lattice: %0.3f deg (%0.3f deg for the line))\n", 
           radians_to_degrees(total_angle), radians_to_degrees(line_angle));
  } else {
    printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(line_angle));
  }

  printf("\nTotal matrix, R, for the line is:\n");
  rmatrix_print(stdout, line_matrix);

  if (args.periodicity > 1) {
    printf("\nTotal matrix, R, for the full system:\n");
    rmatrix_print(stdout, total_matrix);
  }

  if (closed_system) {
    const double R56 = total_matrix[4*BEAM_DOFS + 5];

    printf("\nSynchrotron radiation integrals:\n");
    printf("\tI_1 = %+0.3e  (%+0.3e for the line)\n", I_1, I_1 / args.periodicity);
    printf("\tI_2 = %+0.3e  (%+0.3e for the line)\n", I_2, I_2 / args.periodicity);
    printf("\tI_3 = %+0.3e  (%+0.3e for the line)\n", I_3, I_3 / args.periodicity);
    printf("\tI_4 = %+0.3e  (%+0.3e for the line)\n", I_4, I_4 / args.periodicity);
    printf("\tI_5 = %+0.3e  (%+0.3e for the line)\n", I_5, I_5 / args.periodicity);
    printf("\n");
    printf("x fractional tune:    %f\n", acos(x_trace / 2) / (2*M_PI));
    printf("y fractional tune:    %f\n", acos(y_trace / 2) / (2*M_PI));
    if (gamma_0 > 0)
      printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I_2, gamma_0) / 1e3);
    printf("Momentum compaction:  %0.3e\n", R56 / total_length);
    printf("j_x:                  %+0.3e\n", j_x);
    if (gamma_0 > 0) {
      printf("tau_x:                %0.3f ms\n", 
           1e3 * (2 * args.E_0*1e9 * T_0) / (j_x * e_loss_per_turn(I_2, gamma_0)));
      printf("Natural x emittance:  %0.3f pm.rad\n", 
           1e12 * natural_emittance_x(I_2, I_4, I_5, gamma_0));
      printf("Energy spread:        %0.3e\n", sqrt(energy_spread(I_2, I_3, I_4, gamma_0)));
    }

    arrfree(lin_opt_params.element_etas);
    arrfree(lin_opt_params.element_etaps);
    arrfree(lin_opt_params.element_beta_xs);
    arrfree(lin_opt_params.element_beta_ys);
    arrfree(lin_opt_params.element_curlyH);
    arrfree(lin_opt_params.Ss);
  }

  printf("\n");

  size_t num_turns = 1000 * args.periodicity;
  FILE *aperture_search = fopen("aperture.csv", "w");
  for (double starting_x = -0.03; starting_x < 0.031; starting_x += 0.0005) {
    for (double starting_y = 0.0; starting_y < 0.008; starting_y += 0.0001) {
      double beam[6] = {0};
      beam[0] = starting_x;
      beam[2] = starting_y;
      printf("Tracking %0.3e, %0.3e %zu times", beam[0], beam[2], num_turns);
      size_t turn;
      for (turn=0; turn<num_turns; turn++) {
        track(beam, 1, line, arrlenu(line));
        if ((fabs(beam[0]) > 1.0) || (fabs(beam[2]) > 1.0) || (isnan(beam[0])) || isnan(beam[2])) {
          break;
        }
      }
      printf(": Survived %zu periods\n", turn);
      fprintf(aperture_search, "%0.6e, %0.6e, %zu\n", starting_x, starting_y, turn);
    }
  }
  fclose(aperture_search);

  printf("\n");

  arrfree(line);

  return 0;
}
