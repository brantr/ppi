#include <math.h>
#include <stdlib.h>

void construct_ppi_A(int i, double *q, double *qL, double *qR, double *dq, double *dqVL, double *qh, int n)
{
  double qi   = *q;
  double qimo = *(q-1);
  double qipo = *(q+1);

  double dqi;
  double dqip = qipo - qi;
  double dqim = qi - qimo;

  //eqn 65 MC02
  double dqiVL;

  int ipo = i+1;
  int imo = i-1;

  double ba, bb, bc;
  double osix = 1./6.;


  //first compute qdi, requires q_i-1, q_i+1
  //qdi     = (1/2)( qipo - qimo )
  dq[i]  = 0.5*(q[ipo] - q[imo]); //eqn 64 MC02

  //second, compute qdiVL, requires dq_i, q_i-1, q_i, q_i+1
  //qdiVL   = min(|dqi|, 2|qipo - qi|, 2|qi - qimo|)sgn(dqi)
  dqVL[i] = 0;
  if((q[ipo]-q[i])*(q[i]-q[imo]) > 0)
  {
    ba = 2*fabs(q[ipo]- q[i]);
    bb = 2*fabs(q[i]  - q[imo]);
    dqVL[i] = fmin( fabs(dq[i]), fmin(ba,bb));
    dqVL[i] = copysign(dqVL[i], dq[i]);
  }
}

void construct_ppi_B(int i, double *q, double *qL, double *qR, double *dq, double *dqVL, double *qh, int n)
{
  double qi   = *q;
  double qimo = *(q-1);
  double qipo = *(q+1);

  double dqi;
  double dqip = qipo - qi;
  double dqim = qi - qimo;

  //eqn 65 MC02
  double dqiVL;

  int ipo = i+1;
  int imo = i-1;

  double ba, bb, bc;
  double osix = 1./6.;
  //third, compute q_i+1/2, requires q_i, q_i+1, dq_i+1^VL, dq_i^VL
  //q_i+1/2 = q_i + 0.5*(q_i+1 - q_i) - (1/6)(dq_i+1^VL - dq_i^VL)
  qh[i] = q[i] + 0.5*(q[ipo] - q[i]) - osix*(dqVL[ipo] - dqVL[i]);

  //set q_L,i+1/2 = q_i+1/2, q_R,i+1/2 = q_i+1/2
  qL[i]   = qh[i];  //note qL[i] = q_L,i+1/2
  qR[ipo] = qh[i];  //note qR[i] = q_R,i-1/2
}

void construct_ppi_C(int i, double *q, double *qL, double *qR, double *dq, double *dqVL, double *qh, int n)
{
  double qi   = *q;
  double qimo = *(q-1);
  double qipo = *(q+1);

  double dqi;
  double dqip = qipo - qi;
  double dqim = qi - qimo;

  //eqn 65 MC02
  double dqiVL;

  int ipo = i+1;
  int imo = i-1;

  double ba, bb, bc;
  double osix = 1./6.;

  //double xi_i = 0.0;
  //double xi_i = 0.5;
  //double xi_i = 0.9;
  double xi_i = 1.0;

  double qt;

  //if (q_i+1/2 - q_i)(q_i - q_i-1/2) <=0, then
  //set q_L,i+1/2 = q_i, q_R,i-1/2 = q_i
  ba = qh[i] - q[i];
  bb = q[i]  - qh[imo];
  if(ba*bb<=0)
  {
    qL[i]   = q[i];
    qR[i]   = q[i];
  }

  //if (q_i+1/2 - q_i-1/2)(q_i -0.5(q_i-1/2 +q_i+1/2))>(1/6)(q_i+1/2 - q_i-1/2)^2
  //then q_R,i-1/2 = 3q_i - 2q_i+1/2
  ba = qh[i] - qh[imo];
  bb = q[i]  - 0.5*(qh[imo] + qh[i]);
  bc = qh[i] - qh[imo];
  qt = qR[i];
  if(ba*bb > osix*bc*bc)
    qR[i] = 3*q[i] - 2*qL[i];
    //qR[i] = 3*q[i] - 2*qh[i];


  //if -(1/6)(q_i+1/2 - q_i-1/2)^2 > (q_i+1/2 - q_i-1/2)(q_i - 0.5(q_i+1/2 + q_i-1/2))
  //then q_L,i+1/2 = 3q_i - 2q_i+1/2
  if(-1*osix*bc*bc > ba*bb)
    qL[i] = 3*q[i] - 2*qt;
  //qL[i] = 3*q[i] - 2*qh[i];


  qR[i] = xi_i*qR[i] + (1-xi_i)*q[i];
  qL[i] = xi_i*qL[i] + (1-xi_i)*q[i];

}

double q_ppi(int i, double xi, double *q, double *qL, double *qR)
{
    //return the piecewise parabolic interpolation

  //qi(x) = q_R,i-1/2 + xi(x)[ Delta q_i + q_6i(1 - xi(x))]
  //xi(x) = (x - x_i-1/2)/Delta x_i; x_i-1/2 <= x <= x_i+1/2
  //Delta q_i = q_L,i+1/2 - q_R,i-1/2
  //Delta q_6i = 6[ q_i - (1/2)(q_L,i+1/2 + q_R,i-1/2)]
  double Dq  = qL[i] - qR[i];
  double Dq6 = 6*( q[i] - 0.5*(qL[i] + qR[i]) );
  double qx  = qR[i] + xi*( Dq + Dq6*(1.-xi) );
  return qx;

}