#ifndef NS_ALLOC_H_
#define NS_ALLOC_H_

/**
 * Allocates memory for a rows*cols array of doubles.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 *
 * @param cols Number of columns in the matrix.
 * @param rows Number of rows in the matrix.
 * @return A rows*cols double matrix.
 */
double **alloc_doublematrix(const int cols, const int rows);

/**
 * Allocates memory for a rows*cols array of chars.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 *
 * @param cols Number of columns in the matrix.
 * @param rows Number of rows in the matrix.
 * @return A rows*cols char matrix.
 */
char **alloc_charmatrix(const int cols, const int rows);

/**
 * Frees the memory occupied by a generic matrix.
 *
 * @param m The matrix to be freed.
 */
void free_matrix(void *m);

#endif /* NS_ALLOC_H_ */
