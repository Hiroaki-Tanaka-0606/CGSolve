
void CGSolve(int* N, double** A, double* b, double* threshold, int* iterMax, int* iteration, double* x, double* alpha, double* beta, double** r, bool* verbose);

void ShiftedCGSolve(int* N, double* alpha, double* beta, double** r, int* iteration, double* s, double* x);
