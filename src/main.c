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
#include <math.h>
#include <getopt.h>
#include <errno.h>
#include <unistd.h>

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
	if(argc != 15) {
		fprintf(stderr, "Error: Not enough arguments.\n");
		return -1;
	}

	int rank, nprocesses, initialized = 0;
	double t, delx, dely;
	int i, j, itersor = 0, ifluid = 0, ibound = 0;
	double res;
	double **u = NULL, **v = NULL, **p = NULL, **rhs = NULL, **f = NULL,
		**g = NULL;
	char **flag = NULL;
	int iters = 0;
	unsigned long checker = 0;
	double checker1 = 0.0;
	int ret = 0;

	/* Width of simulated domain. */
	double xlength = atof(argv[1]);
	/* Height of simulated domain. */
	double ylength = atof(argv[2]);
	/* Number of cells horizontally. */
	int imax = atoi(argv[3]);
	/* Number of cells vertically. */
	int jmax = atoi(argv[4]);
	/* Simulation runtime. */
	double t_end = atof(argv[5]);
	/* Duration of each timestep. */
	double del_t = atof(argv[6]);
	/* Safety factor for timestep control. */
	double tau = atof(argv[7]);
	/* Maximum number of iterations in SOR. */
	int itermax = atoi(argv[8]);
	/* Stopping error threshold for SOR. */
	double eps = atof(argv[9]);
	/* Relaxation parameter for SOR. */
	double omega = atof(argv[10]);
	/* Upwind differencing factor in PDE discretisation. */
	double gamma = atof(argv[11]);
	/* Reynolds number. */
	double Re = atof(argv[12]);
	/* Initial X velocity. */
	double ui = atof(argv[13]);
	/* Initial Y velocity. */
	double vi = atof(argv[14]);

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
		fprintf(stderr, "Error: Couldn't allocate memory for matrices.\n");
		ret = -1;
		goto exit;
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

		if(NS_DEBUG_LEVEL) {
			printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
				iters, t+del_t, del_t, itersor, res, ibound);
		}

		update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);
		apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
	}

exit:
	if(u)
		free_matrix(u);
	if(v)
		free_matrix(v);
	if(f)
		free_matrix(f);
	if(g)
		free_matrix(g);
	if(p)
		free_matrix(p);
	if(rhs)
		free_matrix(rhs);
	if(flag)
		free_matrix(flag);

	return ret;
}
