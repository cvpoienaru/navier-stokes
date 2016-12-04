/**
* Modified slightly by D. Orchard (2010) from the classic code from:
*
* Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
* Numerical Simulation in Fluid Dynamics, SIAM, 1998.
* http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html
*/

#ifndef NS_SIMULATION_H_
#define NS_SIMULATION_H_

struct sim_element {
	int row;
	int column;
	double result;
	double north;
	double south;
	double east;
	double west;
};

struct sim_element* alloc_sim_element(void);

void free_sim_element(struct sim_element **element);

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
	const double Re);

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
	const double dely);

int compute_product_sum(
	double **p,
	char **flag,
	const int imax,
	const int jmax,
	const int rank,
	const int size);

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
	const int ifull);

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
	const double dely);

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
	const double tau);

#endif /* NS_SIMULATION_H_ */
