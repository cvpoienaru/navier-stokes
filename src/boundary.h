#ifndef NS_BOUNDARY_H_
#define NS_BOUNDARY_H_

/**
 * Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * matrix edges.
 *
 * @param u Matrix keeping track of the u velocity.
 * @param v Matrix keeping track of the v velocity.
 * @param flag Matrix used to enforce the boundary conditions for u and v
 * matrices.
 * @param imax X axis limit for the matrices.
 * @param jmax Y axis limit for the matrices.
 * @param ui Initial X axis velocity.
 * @param vi Initial Y axis velocity.
 */
void apply_boundary_conditions(
	double **u,
	double **v,
	char **flag,
	const int imax,
	const int jmax,
	const double ui,
	const double vi);

#endif /* NS_BOUNDARY_H_ */
