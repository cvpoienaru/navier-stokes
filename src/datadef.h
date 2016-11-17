#ifndef NS_DATADEF_H_
#define NS_DATADEF_H_

/** This cell is an obstacle/boundary cell. */
#define C_B 0x0000
/** This obstacle cell has a fluid cell to the north. */
#define B_N 0x0001
/** This obstacle cell has a fluid cell to the south. */
#define B_S 0x0002
/** This obstacle cell has a fluid cell to the west. */
#define B_W 0x0004
/** This obstacle cell has a fluid cell to the east. */
#define B_E 0x0008

/** This obstacle cell has a fluid cell to the north west. */
#define B_NW (B_N | B_W)
/** This obstacle cell has a fluid cell to the south west. */
#define B_SW (B_S | B_W)
/** This obstacle cell has a fluid cell to the north east. */
#define B_NE (B_N | B_E)
/** This obstacle cell has a fluid cell to the south east. */
#define B_SE (B_S | B_E)
/** This obstacle cell has fluid cells in all cardinal points. */
#define B_NSEW (B_N | B_S | B_E | B_W)

/** This cell is a fluid cell. */
#define C_F 0x0010

/**
 * Macros for computing Poisson, denoting whether there is an obstacle cell
 * adjacent to some direction.
 */
#define eps_E ((flag[i + 1][j] & C_F) ? 1 : 0)
#define eps_W ((flag[i - 1][j] & C_F) ? 1 : 0)
#define eps_N ((flag[i][j + 1] & C_F) ? 1 : 0)
#define eps_S ((flag[i][j - 1] & C_F) ? 1 : 0)

void write_ppm(
	double** u,
	double** v,
	double** p,
	char** flag,
	const int imax,
	const int jmax,
	const double xlength,
	const double ylength,
	char* outname,
	const int iters,
	const int freq);

unsigned int simplest_checksum_char(char** in, const int imax, const int jmax);
double simplest_checksum(double** in, const int imax, const int jmax);

#endif /* NS_DATADEF_H_ */

