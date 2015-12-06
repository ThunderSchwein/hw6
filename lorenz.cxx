#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;
void f(double *y, double *k1,double *k2,double *k3,const double dx,const double eta);

int main(void)
{
const double dx = 0.1;
const double xmax = 100.0;
//const double eta = 0.4;
//const double psi0 = 1e-5;
double a = 10, b = 28, c = 8/3.0;
double y[2]; y[0] = 1; y[1] = 1; y[2] = 1;

double k1[3];
double k2[3];
double k3[3];

ofstream ex("runge-kutta.dat");
for(int i = 0; i*dx <= xmax; i++)
{

  f(y,k1,k2,k3,dx,eta);
  y[0] += dx*((1.0/6.0)*k1[0]+(2.0/3.0)*k2[0]+(1.0/6.0)*k3[0]);
  y[1] += dx*((1.0/6.0)*k1[1]+(2.0/3.0)*k2[1]+(1.0/6.0)*k3[1]);
  y[2] += dx*((1.0/6.0)*k1[2]+(2.0/3.0)*k2[2]+(1.0/6.0)*k3[2]);
  ex << i*dx << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;
}
ex.close();

exit(0);
}

void f(double *y, double *k1,double *k2,double *k3,const double dx,const double eta)
{

k1[0] = y[1];
k1[1] = (eta - abs(y[0])*abs(y[0]))*y[0];
k1[2]

k2[0] = y[1]+ dx*0.5*k1[1];
k2[1] = (eta - abs(y[0]+dx*0.5*k1[0])*abs(y[0]+dx*0.5*k1[0]))*(y[0]+dx*0.5*k1[0]);
k2[2]

k3[0] = y[1]- dx*k1[1]+ dx*k2[1]*2.0;
k3[1] = (eta - abs(y[0]- dx*k1[0]+ dx*k2[0]*2.0)*abs(y[0]- dx*k1[0]+ dx*k2[0]*2.0))*(y[0]- dx*k1[0]+ dx*k2[0]*2.0);
k3[2] =

return;
}
