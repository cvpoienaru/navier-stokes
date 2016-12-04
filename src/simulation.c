/**
* Modified slightly by D. Orchard (2010) from the classic code from:
*
* Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
* Numerical Simulation in Fluid Dynamics, SIAM, 1998.
* http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html
*/

#include "datadef.h"
#include "alloc.h"
#include "init.h"
#include "simulation.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>

#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))

extern int *ileft, *iright;
extern int nprocs, proc;

struct sim_element* alloc_sim_element(void)
{
	struct sim_element *element = NULL;
	element = (struct sim_element*)malloc(sizeof(struct sim_element));
	if(!element)
		return element;

	element->row = 0;
	element->column = 0;
	element->result = 0.0;
	element->north = 0.0;
	element->south = 0.0;
	element->east = 0.0;
	element->west = 0.0;

	return element;
}

void free_sim_element(struct sim_element **element)
{
	if(!element || !(*element))
		return;

	free(*element);
	*element = NULL;
}

/**
 * Computation of tentative velocity field (f, g).
 */
void compute_tentative_velocity(
	double **u,
	double **v,
	double **f,
	double **g,
	char **flag,
	const int imax,
	const int jmax,
	const double del_t,
	const double delx,
	const double dely,
	const double gamma,
	const double Re)
{
	int i, j;
	double du2dx, duvdy, duvdx, dv2dy, laplu, laplv;

	for (i = 1; i <= imax - 1; ++i) {
		for (j = 1; j <= jmax; ++j) {
			/* Only if both adjacent cells are fluid cells. */
			if ((flag[i][j] & C_F) && (flag[i + 1][j] & C_F)) {
				du2dx = ((u[i][j] + u[i + 1][j]) * (u[i][j] + u[i + 1][j])
					+ gamma * fabs(u[i][j] + u[i + 1][j])
						* (u[i][j] - u[i + 1][j])
					- (u[i - 1][j] + u[i][j]) * (u[i - 1][j] + u[i][j])
					- gamma * fabs(u[i - 1][j] + u[i][j])
						* (u[i - 1][j] - u[i][j]))
					/ (4.0 * delx);
				duvdy = ((v[i][j] + v[i + 1][j]) * (u[i][j] + u[i][j + 1])
					+ gamma * fabs(v[i][j] + v[i + 1][j])
						* (u[i][j] - u[i][j + 1])
					- (v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] + u[i][j])
					- gamma * fabs(v[i][j - 1] + v[i + 1][j - 1])
						*(u[i][j - 1] - u[i][j]))
					/(4.0 * dely);
				laplu = (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j])
						/ delx / delx
					+ (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1])
						/ dely / dely;

				f[i][j] = u[i][j] + del_t * (laplu / Re - du2dx - duvdy);
			} else {
				f[i][j] = u[i][j];
			}
		}
	}

	for (i = 1; i <= imax; ++i) {
		for (j = 1; j <= jmax - 1; ++j) {
			/* Only if both adjacent cells are fluid cells. */
			if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
				duvdx = ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j])
					+ gamma * fabs(u[i][j] + u[i][j + 1])
						* (v[i][j] - v[i + 1][j])
					- (u[i - 1][j] + u[i - 1][j + 1])
						* (v[i - 1][j] + v[i][j])
					- gamma * fabs(u[i - 1][j] + u[i - 1][j + 1])
						* (v[i - 1][j] - v[i][j]))
					/ (4.0 * delx);
				dv2dy = ((v[i][j] + v[i][j + 1]) * (v[i][j] + v[i][j + 1])
					+ gamma * fabs(v[i][j] + v[i][j + 1])
						* (v[i][j] - v[i][j + 1])
					- (v[i][j - 1] + v[i][j]) * (v[i][j - 1] + v[i][j])
					- gamma * fabs(v[i][j - 1] + v[i][j])
						*(v[i][j - 1] - v[i][j]))
					/ (4.0 * dely);

				laplv = (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j])
						/delx/delx
					+(v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1])
						/ dely / dely;

				g[i][j] = v[i][j] + del_t * (laplv / Re - duvdx - dv2dy);
			} else {
				g[i][j] = v[i][j];
			}
		}
	}

	/* F & G at external boundaries. */
	for (j = 1; j <= jmax; ++j) {
		f[0][j] = u[0][j];
		f[imax][j] = u[imax][j];
	}
	for (i = 1; i <= imax; ++i) {
		g[i][0] = v[i][0];
		g[i][jmax] = v[i][jmax];
	}
}

/**
 * Calculate the right hand side of the pressure equation.
 */
void compute_rhs(
	double **f,
	double **g,
	double **rhs,
	char **flag,
	const int imax,
	const int jmax,
	const double del_t,
	const double delx,
	const double dely)
{
	int i, j;

	for (i = 1; i <= imax; ++i) {
		for (j = 1; j <= jmax; ++j) {
			if (flag[i][j] & C_F) {
				/* Only for fluid and non-surface cells. */
				rhs[i][j] = ((f[i][j] - f[i - 1][j]) / delx
					+ (g[i][j] - g[i][j - 1]) / dely) / del_t;
			}
		}
	}
}

int compute_product_sum(
	double **p,
	char **flag,
	const int imax,
	const int jmax,
	const int rank,
	const int size)
{
	int i, j;
	double p0 = 0.0;
	int sum_chunk_size = (imax * jmax) / size;
	int sum_start = rank * sum_chunk_size;
	int sum_end = (rank == size - 1)
		? imax * jmax
		: (rank + 1) * sum_chunk_size;
	int sum_start_j = (sum_start % jmax) + 1;
	int sum_start_i = ((sum_start - sum_start_j) % imax) + 1;
	int sum_end_j = (sum_end % jmax) + 1;
	int sum_end_i = ((sum_end - sum_end_j) % imax) + 1;
	MPI_Status stat;

	/* Calculate sum of squares */
	i = sum_start_i;
	j = sum_start_j;
	while(i <= imax) {
		while(j <= jmax) {
			if(i == sum_end_i && j == sum_end_j)
				goto mpi_send;

			if(flag[i][j] & C_F)
				p0 += p[i][j] * p[i][j];
			++j;
		}
		j = 1;
		++i;
	}

mpi_send:
	if(rank == MASTER) {
		int recv = size - 1;
		double psum = 0.0;

		while(recv > 0) {
			MPI_Recv(&psum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
				MPI_COMM_WORLD, &stat);
			--recv;
			p0 += psum;
		}
		for(i = 1; i < size; ++i) {
			MPI_Send(&p0, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
		}
	} else {
		MPI_Send(&p0, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
		MPI_Recv(&p0, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
			MPI_COMM_WORLD, &stat);
	}

	return p0;
}

/**
 * Red/Black SOR to solve the poisson equation.
 */
int poisson(
	double **p,
	double **rhs,
	char **flag,
	const int imax,
	const int jmax,
	const double delx,
	const double dely,
	const double eps,
	const int itermax,
	const double omega,
	double *res,
	const int ifull)
{
	int rank, size, rb;
	int i, j, iter;
	double add, beta_2, beta_mod;
	double p0 = 0.0;
	double rdx2 = 1.0 / (delx * delx);
	double rdy2 = 1.0 / (dely * dely);
	beta_2 = -omega / (2.0 * (rdx2 + rdy2));

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	p0 = compute_product_sum(p, flag, imax, jmax, rank, size);
	MPI_Barrier(MPI_COMM_WORLD);

	p0 = sqrt(p0 / ifull);
	if (p0 < 0.0001)
		p0 = 1.0;

	/* Red/Black SOR-iteration. */
	for (iter = 0; iter < itermax; ++iter) {
		for (rb = 0; rb <= 1; ++rb) {
			for (i = 1; i <= imax; ++i) {
				for (j = 1; j <= jmax; ++j) {
					if ((i + j) % 2 != rb)
						continue;

					if (flag[i][j] == (C_F | B_NSEW)) {
						/* Five point star for interior fluid cells. */
						p[i][j] = (1.0 - omega) * p[i][j] - beta_2
							* ((p[i + 1][j] + p[i - 1][j]) * rdx2
								+ (p[i][j + 1] + p[i][j - 1]) * rdy2
								- rhs[i][j]);
					} else if (flag[i][j] & C_F) {
						/* Modified star near boundary. */
						beta_mod = -omega / ((eps_E + eps_W) * rdx2
							+(eps_N + eps_S) * rdy2);
						p[i][j] = (1.0 - omega) * p[i][j] - beta_mod
							* ((eps_E * p[i + 1][j] + eps_W * p[i - 1][j])
									* rdx2
								+ (eps_N * p[i][j + 1] + eps_S * p[i][j - 1])
									* rdy2
							- rhs[i][j]);
					}
				}
			}
		}

		/* Partial computation of residual. */
		*res = 0.0;
		for (i = 1; i <= imax; ++i) {
			for (j = 1; j <= jmax; ++j) {
				if (flag[i][j] & C_F) {
					/* Only fluid cells. */
					add = (eps_E * (p[i + 1][j] - p[i][j])
						- eps_W * (p[i][j] - p[i - 1][j])) * rdx2
						+ (eps_N * (p[i][j + 1] - p[i][j])
							- eps_S * (p[i][j] - p[i][j - 1])) * rdy2
						- rhs[i][j];
					*res += add * add;
				}
			}
		}
		*res = sqrt((*res) / ifull) / p0;

		/* Convergence ? */
		if (*res < eps)
			break;
	}

	return iter;
}

/**
 * Update the velocity values based on the tentative velocity values and the new
 * pressure matrix.
*/
void update_velocity(
	double **u,
	double **v,
	double **f,
	double **g,
	double **p,
	char **flag,
	const int imax,
	const int jmax,
	const double del_t,
	const double delx,
	const double dely)
{
	int i, j;

	for (i = 1; i <= imax - 1; ++i) {
		for (j = 1; j <= jmax; ++j) {
			/* Only if both adjacent cells are fluid cells. */
			if ((flag[i][j] & C_F) && (flag[i + 1][j] & C_F)) {
				u[i][j] = f[i][j] - (p[i + 1][j] - p[i][j]) * del_t / delx;
			}
		}
	}

	for (i = 1; i <= imax; ++i) {
		for (j = 1; j <= jmax - 1; ++j) {
			/* Only if both adjacent cells are fluid cells. */
			if ((flag[i][j] & C_F) && (flag[i][j + 1] & C_F)) {
				v[i][j] = g[i][j] - (p[i][j + 1] - p[i][j]) * del_t / dely;
			}
		}
	}
}

/**
 * Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
void set_timestep_interval(
	double *del_t,
	const int imax,
	const int jmax,
	const double delx,
	const double dely,
	double **u,
	double **v,
	const double Re,
	const double tau)
{
	int i, j;
	double umax, vmax, deltu, deltv, deltRe;

	/* del_t satisfying CFL conditions. */
	if (tau >= 1.0e-10) {
		/* Else no time stepsize control. */
		umax = 1.0e-10;
		vmax = 1.0e-10;
		for (i = 0; i <= imax + 1; ++i) {
			for (j = 1; j <= jmax + 1; ++j) {
				umax = max(fabs(u[i][j]), umax);
			}
		}

		for (i = 1; i <= imax + 1; ++i) {
			for (j = 0; j <= jmax + 1; ++j) {
				vmax = max(fabs(v[i][j]), vmax);
			}
		}

		deltu = delx / umax;
		deltv = dely / vmax; 
		deltRe = 1 / (1 / (delx * delx) + 1 / (dely * dely)) * Re / 2.0;

		if (deltu < deltv) {
			*del_t = min(deltu, deltRe);
		} else {
			*del_t = min(deltv, deltRe);
		}

		/* Multiply by safety factor. */
		*del_t = tau * (*del_t);
	}
}
