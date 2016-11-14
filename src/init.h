#ifndef NS_INIT_H_
#define NS_INIT_H_

void load_flag_from_pgm(char **flag, int imax, int jmax, char *filename);
void init_flag(char **flag, int imax, int jmax, double delx, double dely,
    int *ibound);

#endif /* NS_INIT_H_ */

