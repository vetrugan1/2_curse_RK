#include <iostream>
#include <math.h>

using namespace std;

void right(double *x, double *f, double t) {
	f[0] = (double)x[1];
	f[1] = (double)0.5 * x[0] + (double)0.5 * x[1] + 3 * exp(t / 2);
}

void Eiler(double t, double h, double* x) {
	double xr[2], f[2];
	int i;

	right(x, f, t);

	for (i = 0; i < 2; i++) {
		xr[i] = (double)x[i] + h * f[i];
	}
	for (i = 0; i < 2; i++) {
		x[i] = (double)x[i] + h * f[i];
	}
}

// Method Runge-Kutta-2
void RK_2(double t, double h, double* x) {
	double xr[2], f[2];
	int i;
	double h2;

	h2 = (double)0.5 * h;

	right(x, f, t);
	for (i = 0; i < 2; i++) {
		xr[i] = (double)x[i] + h2 * f[i];
	}

	right(xr, f, t + h2);
	for (i = 0; i < 2; i++) {
		x[i] = (double)x[i] + h * f[i];
	}
}

// Method Runge-Kutta-4
void RK_4(double t, double h, double* x) {
	double xr[2], f1[2], f2[2], f3[2], f4[2];
	int i;
	double h2, h6;

	h2 = (double)0.5 * h;
	h6 = (double)h / 6;

	// first 
	right(x, f1, t);
	for (i = 0; i < 2; i++) {
		xr[i] = (double)x[i] + h2 * f1[i];
	}

	//second
	right(xr, f2, t+h2);
	for (i = 0; i < 2; i++) {
		xr[i] = (double)x[i] + h2 * f2[i];
	}

	//third
	right(xr, f3, t + h2);
	for (i = 0; i < 2; i++) {
		xr[i] = (double)x[i] + h * f3[i];
	}

	//fourth
	right(xr, f4, t + h2);
	for (i = 0; i < 2; i++) {
		x[i] = (double)x[i] + h6 * (f1[i]+2*(f2[i]+f3[i])+f4[i]);
	}
}

int main()
{
	double xk1[2], xk2[2], xk3[2];
	double t = 0, d2h = 1;
	double step = 0.01;
	double tr = 0;

	xk1[0] = -4;
	xk2[0] = -4;
	xk3[0] = -4;

	xk1[1] = -2.5;
	xk2[1] = -2.5;
	xk3[1] = -2.5;


	while (t < d2h) {
		tr = exp(t) - 6 * exp(t / 2) + exp(-t / 2);
		Eiler(t, step, xk1);
		RK_2(t, step, xk2);
		RK_4(t, step, xk3);
		t += step;
		cout << "xk1[0] = " << xk1[0] << "\t"
			<< "xk2[0] = " << xk2[0] << "\t" 
			<< "xk3[0] = " << xk3[0] << "\t" 
			<<  "tr = " << tr << endl;
	}
}
