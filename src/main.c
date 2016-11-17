/**
* Modified slightly by D. Orchard (2010) from the classic code from:
*
* Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
* Numerical Simulation in Fluid Dynamics, SIAM, 1998.
* http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html
*/

#include "datadef.h"
#include "alloc.h"
#include "boundary.h"
#include "init.h"
#include "simulation.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <errno.h>

/* Functions used for comparing computations when debugging. */
unsigned int simplest_checksum_char(char** in, int imax, int jmax)
{
	int i;
	int j;
	unsigned int checksum = 0;

	for (i = 0; i < imax + 2; ++i) {
		for (j = 0; j < jmax + 2; ++j) {
			checksum += in[i][j] * i;
		}
	}

	return checksum;
}

double simplest_checksum(double** in, int imax, int jmax)
{
	int i;
	int j;
	double checksum = 0.0;

	for (i = 0; i < imax + 2; ++i) {
		for (j = 0; j < jmax + 2; ++j) {
			checksum += in[i][j] * ((double)(i * jmax) + j);
		}
	}

	return checksum;
}

int main(int argc, char **argv)
{
	/* Verbosity level. */
	int verbose = 1;
	/* Width of simulated domain. */
	double xlength = 22.0;
	/* Height of simulated domain. */
	double ylength = 4.1;
	/* Number of cells horizontally. */
	int imax = 660;
	/* Number of cells vertically. */
	int jmax = 120;

	char *outname;
	int output = 0;
	int output_frequency = 0;

	/* Simulation runtime. */
	double t_end = 40;
	/* Duration of each timestep. */
	double del_t = 0.3;
	/* Safety factor for timestep control. */
	double tau = 0.5;

	/* Maximum number of iterations in SOR. */
	int itermax = 100;
	/* Stopping error threshold for SOR. */
	double eps = 0.001;
	/* Relaxation parameter for SOR. */
	double omega = 1.7;
	/* Upwind differencing factor in PDE discretisation. */
	double gamma = 0.9;

	/* Reynolds number. */
	double Re = 150.0;
	/* Initial X velocity. */
	double ui = 1.0;
	/* Initial Y velocity. */
	double vi = 0.0;

	double t, delx, dely;
	int i, j, itersor = 0, ifluid = 0, ibound = 0;
	double res;
	double **u, **v, **p, **rhs, **f, **g;
	char **flag;
	int init_case, iters = 0;

	unsigned long checker = 0;
	double checker1 = 0.0;

	if (argc > 1) {
		output = 1;
		outname = argv[1];
		output_frequency = 1;
	}

	if (argc > 2) {
		output_frequency = atoi(argv[2]);
	}

	delx = xlength/imax;
	dely = ylength/jmax;

	/* Allocate arrays. */
	u = alloc_doublematrix(imax + 2, jmax + 2);
	v = alloc_doublematrix(imax + 2, jmax + 2);
	f = alloc_doublematrix(imax + 2, jmax + 2);
	g = alloc_doublematrix(imax + 2, jmax + 2);
	p = alloc_doublematrix(imax + 2, jmax + 2);
	rhs = alloc_doublematrix(imax + 2, jmax + 2);
	flag = alloc_charmatrix(imax + 2, jmax + 2);

	if (!u || !v || !f || !g || !p || !rhs || !flag) {
		fprintf(stderr, "Couldn't allocate memory for matrices.\n");
		return 1;
	}

	/* Set up initial values. */
	for (i = 0; i <= imax + 1; ++i) {
		for (j = 0; j <= jmax + 1; ++j) {
			checker += (i * jmax) + j + 1;
			checker1 += (i * jmax) + j + 1.0;
			u[i][j] = ui;
			v[i][j] = vi;
			p[i][j] = 0.0;
		}
	}

	init_flag(flag, imax, jmax, delx, dely, &ibound);
	apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);

	for (t = 0.0; t < t_end; t += del_t, ++iters) {
		set_timestep_interval(&del_t, imax, jmax, delx, dely, u, v, Re, tau);

		ifluid = (imax * jmax) - ibound;
		compute_tentative_velocity(u, v, f, g, flag, imax, jmax,
			del_t, delx, dely, gamma, Re);
		compute_rhs(f, g, rhs, flag, imax, jmax, del_t, delx, dely);

		itersor = 0;
		if (ifluid > 0) {
			itersor = poisson(p, rhs, flag, imax, jmax, delx, dely,
				eps, itermax, omega, &res, ifluid);
		}

		/* printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
		iters, t+del_t, del_t, itersor, res, ibound); */

		update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);
		apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);

		if (output && (iters % output_frequency == 0)) {
			write_ppm(u, v, p, flag, imax, jmax, xlength, ylength, outname,
				iters, output_frequency);
		}
	}

	free_matrix(u);
	free_matrix(v);
	free_matrix(f);
	free_matrix(g);
	free_matrix(p);
	free_matrix(rhs);
	free_matrix(flag);

	return 0;
}
