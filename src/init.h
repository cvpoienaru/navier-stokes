#ifndef NS_INIT_H_
#define NS_INIT_H_

/**
 * Initializes the flag array, marking any obstacle cells and edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 *
 * @param flag The flag matrix used to enforce boundary conditions.
 * @param imax X axis limit for the flag matrix.
 * @param jmax Y axis limit for the flag matrix.
 * @param delx Cell X axis size.
 * @param dely Cell Y axis size.
 * @param ibound The ibound.
 */
void init_flag(
	char **flag,
	const int imax,
	const int jmax,
	const double delx,
	const double dely,
	int *ibound);

// DO WE NEED THIS ??
void load_flag_from_pgm(char **flag, int imax, int jmax, char *filename);

#endif /* NS_INIT_H_ */

