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
double a = 10, b = 28, c = 8/3.0;
double y[2]; y[0] = 1; y[1] = 1; y[2] = 1;

double k1[3];
double k2[3];
double k3[3];

ofstream Ausgabe("runge-kutta.dat");
for(int i = 0; i*dx <= xmax; i++)
{
  f(y,k1,k2,k3,k4,dx,eta);
  y[0] += dx*((1.0/6.0)*k1[0]+(2.0/3.0)*k2[0]+(2.0/6.0)*k3[0]+(1.0/6.0)*k4[0]);
  y[1] += dx*((1.0/6.0)*k1[1]+(2.0/3.0)*k2[1]+(2.0/6.0)*k3[1]+(1.0/6.0)*k4[1]);
  y[2] += dx*((1.0/6.0)*k1[2]+(2.0/3.0)*k2[2]+(2.0/6.0)*k3[2]+(1.0/6.0)*k4[2]);
  Ausgabe << i*dx << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;
}
Ausgabe.close();

exit(0);
}

void f(double *y, double *k1,double *k2,double *k3, double *k4, const double dx,const double eta)
{

k1[0] = a*(y[1]-y[0]); //a*(y-x)
k1[1] = y[0]*(b-y[2]) - y[1]; //x*(b-z) - y
k1[2] = y[0]*y[1] - c*y[2]; //x*y-c*z

k2[0] = (y[1]+k1[1]*dt/2 -(y[0]+k1[0]*dt/2))*a;
k2[1] = (y[0]+k1[0]*dt/2)*(b-(y[2]+k1[2]*dt/2)) - (y[1] + k1[1]*dt/2); 
k2[2] = (y[0]+k1[0]*dt/2)*(y[1] + k1[1]*dt/2) - c*(y[2]+k1[2]*dt/2);

k3[0] = (y[1]+k2[1]*dt/2  -(y[0]+k2[0]*dt/2))*a;
k3[1] = (y[0]+k2[0]*dt/2) *(b-(y[2]+k2[2]*dt/2)) - (y[1] + k2[1]*dt/2); 
k3[2] = (y[0]+k2[0]*dt/2) *(y[1] + k2[1]*dt/2) - c*(y[2]+k2[2]*dt/2);

k4[0] = (y[1]+k3[1]*dt  -(y[0]+k3[0]*dt))*a;
k4[1] = (y[0]+k3[0]*dt) *(b-(y[2]+k3[2]*dt)) - (y[1] + k3[1]*dt); 
k4[2] = (y[0]+k3[0]*dt) *(y[1] + k3[1]*dt) - c*(y[2]+k3[2]*dt);

return;
}
