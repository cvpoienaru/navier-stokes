/**
* Modified slightly by D. Orchard (2010) from the classic code from:
*
* Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
* Numerical Simulation in Fluid Dynamics, SIAM, 1998.
* http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "datadef.h"
#include "alloc.h"
#include "init.h"
#include "simulation.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifdef __cplusplus
}
#endif

#include "tbb/tbb.h"
#include "tbb/parallel_for.h"

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
	int i, j, iter, rank, size;
	double add, beta_2;
	double p0 = 0.0;
	double rdx2 = 1.0 / (delx * delx);
	double rdy2 = 1.0 / (dely * dely);
	beta_2 = -omega / (2.0 * (rdx2 + rdy2));

	/* Red-black value. */
	int rb;

	/* Calculate sum of squares */
	for (i = 1; i <= imax; ++i) {
		for (j = 1; j <= jmax; ++j) {
			if (flag[i][j] & C_F)
				p0 += p[i][j] * p[i][j];
		}
	};

	p0 = sqrt(p0 / ifull);
	if (p0 < 0.0001)
		p0 = 1.0;

	int nthread = 32;
	tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);

	/* Red/Black SOR-iteration. */
	for (iter = 0; iter < itermax; ++iter) {
		for (rb = 0; rb <= 1; ++rb) {
			tbb::parallel_for(1, imax, [&](int k) {
				int l;
				for(l = 1; l <= jmax; ++l) {
					if ((k + l) % 2 == rb) {
						if (flag[k][l] == (C_F | B_NSEW)) {
							p[k][l] = (1.0 - omega) * p[k][l] - beta_2
								* ((p[k + 1][l] + p[k - 1][l]) * rdx2
									+ (p[k][l + 1] + p[k][l - 1]) * rdy2
									- rhs[k][l]);
						} else if (flag[k][l] & C_F) { 
							double beta_mod = -omega / ((eps_E(k, l) + eps_W(k, l)) * rdx2
								+(eps_N(k, l) + eps_S(k, l)) * rdy2);
							p[k][l] = (1.0 - omega) * p[k][l] - beta_mod
								* ((eps_E(k, l) * p[k + 1][l] + eps_W(k, l) * p[k - 1][l])
										* rdx2
									+ (eps_S(k, l) * p[k][l + 1] + eps_N(k, l) * p[k][l - 1])
										* rdy2
								- rhs[k][l]);
						}
					}
				}
			});
		}

		/* Partial computation of residual. */
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
		};
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
