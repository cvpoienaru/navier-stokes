#ifndef NS_BOUNDARY_H_
#define NS_BOUNDARY_H_

void apply_boundary_conditions(double **u, double **v, char **flag,
    const int imax, const int jmax, const double ui, const double vi);

#endif /* BOUNDARY_H_ */

