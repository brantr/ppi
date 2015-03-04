#ifndef PPI_H
#define PPI_H
void construct_ppi_A(int i, double *q, double *qL, double *qR,double *dq, double *dqVL, double *qh, int n); 
void construct_ppi_B(int i, double *q, double *qL, double *qR,double *dq, double *dqVL, double *qh, int n); 
void construct_ppi_C(int i, double *q, double *qL, double *qR,double *dq, double *dqVL, double *qh, int n); 
double q_ppi(int ix, double xi, double *q, double *qL, double *qR);
#endif //PPI_H