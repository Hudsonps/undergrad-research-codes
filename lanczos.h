#ifndef _LANCZOS_
#define _LANCZOS_
const int NBAND = 100;
const double TOLERANCE = 1.e-12;
const double THRESHOLD = 1.e-13;
const double E_MIN = 1.e-7;
double ksi(int n, int dim = NBAND);
#endif
