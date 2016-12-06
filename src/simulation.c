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

static const int compute_product_sum(
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

static void range_to_index(
	const int cols,
	const int s_i,
	const int s_j,
	const int e_n,
	int *e_i,
	int *e_j)
{
	int n = e_n + (s_i * cols + s_j);
	*e_j = n % cols;
	*e_i = n / cols;
}

static void index_to_range(
	const int cols,
	const int s_i,
	const int s_j,
	const int e_i,
	const int e_j,
	int *e_n)
{
	*e_n = (e_i * cols + e_j) - (s_i * cols + s_j);
}

static double* copy_poisson_range(
	double **p,
	const int rows,
	const int cols,
	const int start_i,
	const int start_j,
	const int end_i,
	const int end_j,
	int *elements_no)
{
	int i = start_i;
	int j = start_j;
	int k = -1;
	double *range = NULL;

	*elements_no = ((end_i - start_i) * cols + 1) + (end_j - start_j);
	range = (double*)malloc(*elements_no * sizeof(double));
	if(!range)
		goto exit;

	while(i < rows) {
		while(j < cols) {
			if(i == rows - 1 && j == cols)
				goto exit;

			range[++k] = p[i][j];
			++j;
		}

		j = 0;
		++i;
	}

exit:
	return range;
}

static void do_slave_p(
	double **rhs,
	char **flag,
	const int rows,
	const int cols,
	const double omega,
	const double beta_2,
	const double rdx2,
	const double rdy2)
{
	int i, valid, range_size, rb;
	int e_i, e_j, s_i, s_j;
	int north_i, south_i, east_i, west_i;
	double north, south, east, west;
	double beta_mod;
	double *range = NULL;
	MPI_Status stat;

	MPI_Recv(&s_i, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD,
		&stat);
	MPI_Recv(&s_j, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD,
		&stat);
	MPI_Recv(&range_size, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD,
		&stat);
	MPI_Recv(&rb, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD,
		&stat);

	range = (double*)malloc(range_size * sizeof(double));
	if(!range)
		return;

	MPI_Recv(range, range_size, MPI_DOUBLE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD,
		&stat);

	for(i = cols; i < range_size - cols; ++i) {
		range_to_index(cols, s_i, s_j, i, &e_i, &e_j);

		valid = (e_i >= 1 && e_i <= rows - 2) && (e_j >= 1 && e_j <= cols - 2)
			&& ((e_i + e_j) % 2 == rb);
		if(!valid)
			continue;

		index_to_range(cols, s_i, s_j, e_i - 1, e_j, &north_i);
		index_to_range(cols, s_i, s_j, e_i + 1, e_j, &south_i);
		index_to_range(cols, s_i, s_j, e_i, e_j + 1, &east_i);
		index_to_range(cols, s_i, s_j, e_i, e_j - 1, &west_i);
		north = range[north_i];
		south = range[south_i];
		east = range[east_i];
		west = range[west_i];

		if (flag[e_i][e_j] == (C_F | B_NSEW)) {
			/* Five point star for interior fluid cells. */
			range[i] = (1.0 - omega) * range[i] - beta_2
				* ((south + north) * rdx2
					+ (east + west) * rdy2
					- rhs[e_i][e_j]);
		} else if (flag[e_i][e_j] & C_F) {
			/* Modified star near boundary. */
			beta_mod = -omega / ((eps_E(e_i, e_j) + eps_W(e_i, e_j)) * rdx2
				+ (eps_N(e_i, e_j) + eps_S(e_i, e_j)) * rdy2);
			range[i] = (1.0 - omega) * range[i] - beta_mod
				* ((eps_E(e_i, e_j) * south + eps_W(e_i, e_j) * north) * rdx2
					+ (eps_S(e_i, e_j) * east + eps_N(e_i, e_j) * west) * rdy2
					- rhs[e_i][e_j]);
		}
	}

	MPI_Send(&range_size, 1, MPI_INT, MASTER, MASTER, MPI_COMM_WORLD);
	MPI_Send(range, range_size, MPI_DOUBLE, MASTER, MASTER, MPI_COMM_WORLD);
	free(range);
}

static void do_master_p(
	double **p,
	const int rows,
	const int cols,
	int rb,
	const int size)
{
	int curr_even_chunk, curr_odd_chunk, l, source;
	int i, j, k, proc;
	int n = rows;
	int m = cols;
	int evens = ((n * m) - 2 * (n + m) + 4) / 2 + (n * m % 2 == 0 ? 0 : 1);
	int odds = ((n * m) - 2 * (n + m) + 4) / 2;
	int even_chunk_size = evens / (size - 1);
	int odd_chunk_size = odds / (size - 1);
	int s_i, s_j, e_i, e_j, first = TRUE;
	int range_size = 0;
	int *proc_range_size = NULL;
	double *range = NULL;
	double **proc_range = NULL;
	MPI_Status stat;

	proc = 0;
	k = 0;
	for(i = 1; i <= rows - 2; ++i) {
		for(j = 1; j <= cols - 2; ++j) {
			curr_even_chunk = (proc == size - 2)
				? evens - proc * even_chunk_size
				: even_chunk_size;
			curr_odd_chunk = (proc == size - 2)
				? odds - proc * odd_chunk_size
				: odd_chunk_size;

			if((i + j) % 2 != rb)
				continue;

			if(first == TRUE) {
				s_i = i - 1;
				s_j = j;
				first = FALSE;
			}

			++k;
			if(k == (rb == 0 ? curr_even_chunk : curr_odd_chunk)) {
				e_i = i + 1;
				e_j = j;
				k = 0;
				first = TRUE;
				++proc;

				range = copy_poisson_range(p, n, m, s_i, s_j, e_i, e_j,
					&range_size);

				MPI_Send(&s_i, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
				MPI_Send(&s_j, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
				MPI_Send(&range_size, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
				MPI_Send(&rb, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
				MPI_Send(range, range_size, MPI_DOUBLE, proc, proc,
					MPI_COMM_WORLD);
				free(range);
			}
		}
	}

	proc_range = (double**)malloc((size - 1) * sizeof(double*));
	if(!proc_range)
		return;

	proc_range_size = (int*)malloc((size - 1) * sizeof(int));
	if(!proc_range_size)
		return;

	proc = size - 1;
	while(proc > 0) {
		MPI_Recv(&range_size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
			MPI_COMM_WORLD, &stat);

		source = stat.MPI_SOURCE;
		proc_range[source - 1] = (double*)malloc(range_size * sizeof(double));
		proc_range_size[source - 1] = range_size;
		if(!proc_range[source - 1])
			return;

		MPI_Recv(proc_range[source - 1], range_size, MPI_DOUBLE, source,
			MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
		--proc;
	}

	l = 0;
	k = cols;
	for(i = 1; i <= rows - 2; ++i) {
		for(j = 1; j <= cols - 2; ++j) {
			if((i + j) % 2 != rb)
				continue;

			/*while(k < ) {
				//if(range[k])
			}*/
		}
	}

	/*k = -1;
	proc = 0;
	first = TRUE;
	for(i = 1; i <= rows - 2; ++i) {
		for(j = 1; j <= cols - 2; ++j) {
			if(k == proc_range_size[proc]) {
				k = -1;
				first = TRUE;
				++proc;
			}
			if(k != -1)
				++k;

			if((i + j) % 2 != rb)
				continue;

			if(first == TRUE) {
				first = FALSE;
				k = cols;
			}

			p[i][j] = proc_range[proc][k];
		}
	}*/

	for(i = 0; i < size - 1; ++i)
		free(proc_range[i]);
	free(proc_range);
	free(proc_range_size);
}

static void compute_p(
	double **p,
	double **rhs,
	char **flag,
	const int imax,
	const int jmax,
	const double omega,
	const double beta_2,
	const double rdx2,
	const double rdy2,
	const int rank,
	const int size)
{
	int rb;
	int n = imax + 2;
	int m = jmax + 2;

	for(rb = 0; rb <= 1; ++rb) {
		(rank == MASTER)
			? do_master_p(p, n, m, rb, size)
			: do_slave_p(rhs, flag, n, m, omega, beta_2, rdx2, rdy2);
	}
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
	int rank, size;
	int i, j, iter;
	double add, beta_2;
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
		compute_p(p, rhs, flag, imax, jmax, omega, beta_2, rdx2, rdy2, rank,
			size);
		MPI_Barrier(MPI_COMM_WORLD);

		// Partial computation of residual.
		*res = 0.0;
		for (i = 1; i <= imax; ++i) {
			for (j = 1; j <= jmax; ++j) {
				if (flag[i][j] & C_F) {
					/* Only fluid cells. */
					add = (eps_E(i, j) * (p[i + 1][j] - p[i][j])
						- eps_W(i, j) * (p[i][j] - p[i - 1][j])) * rdx2
						+ (eps_S(i, j) * (p[i][j + 1] - p[i][j])
							- eps_N(i, j) * (p[i][j] - p[i][j - 1])) * rdy2
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
