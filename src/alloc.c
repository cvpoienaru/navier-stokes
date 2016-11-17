#include "alloc.h"

#include <stdlib.h>

/**
 * Allocates memory for a rows*cols array of doubles.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 *
 * @param cols Number of columns in the matrix.
 * @param rows Number of rows in the matrix.
 * @return A rows*cols double matrix.
 */
double** alloc_doublematrix(const int cols, const int rows)
{
	int i;
	double **m;

	m = (double**)malloc(cols * sizeof(double*));
	if (!m)
		return NULL;

	double *els = (double*)calloc(rows * cols, sizeof(double));
	if (!els)
		return NULL;

	for (i = 0; i < cols; ++i) {
		m[i] = &els[rows * i];
	}

	return m;
}

/**
 * Allocates memory for a rows*cols array of chars.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 *
 * @param cols Number of columns in the matrix.
 * @param rows Number of rows in the matrix.
 * @return A rows*cols char matrix.
 */
char** alloc_charmatrix(const int cols, const int rows)
{
	int i;
	char **m;

	m = (char**)malloc(cols * sizeof(char*));
	if (!m)
		return NULL;

	char *els = (char*)malloc(rows * cols * sizeof(char));
	if (!els)
		return NULL;

	for (i = 0; i < cols; ++i) {
		m[i] = &els[rows * i];
	}

	return m;
}

/**
 * Frees the memory occupied by a generic matrix.
 *
 * @param m The matrix to be freed.
 */
void free_matrix(void *m)
{
	void **els = (void**)m;
	free(els[0]);	/* Deallocate the block of array elements */
	free(m);		/* Deallocate the block of column pointers */
}
