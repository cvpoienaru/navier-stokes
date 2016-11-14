#ifndef NS_ALLOC_H_
#define NS_ALLOC_H_

double **alloc_doublematrix(const int cols, const int rows);
char **alloc_charmatrix(const int cols, const int rows);
void free_matrix(void *m);

#endif /* NS_ALLOC_H_ */
