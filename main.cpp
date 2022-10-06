#include <iostream>
#include <math.h>
#include <cstring>
//#include "data.h"

using namespace std;

//=====================
const double pi = 3.141592653589793;

int N = 20;

int flag = 0;
double err = 1;
double eps = 0.05;

double *prevd = new double[N + 1];
double *nextd = new double[N + 1];
double *prehalfd = new double[N + 1];
double *nehalfd = new double[N + 1];

double space_border = 1;
double h_step = space_border / N; 

int T = 1.1 * 2 * N * N;
double time_border = 1;
double tau = time_border / T;


//=====================

void update() {
	N = N * 2;
	
	flag = 0;
	
	delete[] prevd;
	delete[] nextd;
	delete[] prehalfd;
	delete[] nehalfd;
	
	prevd = new double[N + 1];
	nextd = new double[N + 1];
	prehalfd = new double[N + 1];
	nehalfd = new double[N + 1];

	h_step = space_border / N; 

	T = 1.1 * 2 * N * N;
	tau = time_border / T;
	}

double solution(double x, double t) {
	return exp(-pi*pi*t)*sin(pi*x);
}

double U_0 (double x) { /* U(x,0) = U_0(x)*/
	return sin(pi*x);
}

double U_left(double t) { /* U(0,t) = U_left(t)*/
	return 0;
}

double U_right(double t) { /* U(l,t) = U_right(t)*/ 
	return 0;
}

double f(double x, double t) {
	return 0;
}

/* y{i,j+1} = y{i,j} + tau/h^2 *(y{i+1,j} - 2y{i,j} + y{i-1,j}) + tau * fi */
int explicit_scheme() {
	cout << "N=" << N << endl;
	int k = 0;
	double max = 0;
	
	FILE *output;
	output = fopen("ex_exact.txt", "w");
	FILE *output_1;
	output_1 = fopen("ex_solution.txt", "w");

	for (int i = 0; i < N + 1; i++) {
		prevd[i] = U_0(i * h_step);
		prehalfd[i] = U_0(h_step / 2 + i * h_step);
	}

	for (int j = 1; j < T + 1 ; j++) {
		//cout << "timestep = " << j << endl;
		for (int i = 1; i < N; i++) {
			nextd[i] = prevd[i] + (tau / (h_step * h_step))*(prevd[i + 1] - 2 * prevd[i] + prevd[i - 1]);
			nehalfd[i] = prehalfd[i] + (tau / (h_step * h_step))*(prehalfd[i + 1] - 2 * prehalfd[i] + prehalfd[i - 1]);
		}

		nextd[0] = 0;
		nextd[N] = 0;

		for (int i = 0; i < N + 1; i++)  {
			prevd[i] = nextd[i];
			prehalfd[i] = nehalfd[i];
		}
		
		for (int i = 0; i < N + 1; i++) {
			if (max < abs(prevd[i] - prehalfd[i]) * 2) {
				max = abs(prevd[i] - prehalfd[i]) * 2;
			}	
		}
		
		if (max > eps) {
			cout << max << "\nLets update the grid" << endl;
			update();
			fclose(output);
			fclose(output_1);
			return 0;
		}

		if (k%500 == 0) {
			for (int i = 0; i < N + 1; i++) {
				fprintf(output, "%f %f\n", i*h_step, solution(i*h_step, j*tau));
				fprintf(output_1, "%f %f\n", i*h_step, prevd[i]);	
				//cout << "step " << i*h_step << " our_solution=" << prevd[i] << " solution " << solution(i*h_step, j*tau) << " delta " << abs(solution(i*h_step, j*tau)-prevd[i]) << endl;
			}
		}
				
		k++;
	}
	cout << "max=" << max << endl;
	return 1;
}

void progonka(double time_step, double* preArray, double* neArray) {
	double *alpha = new double[N + 1];
	double *beta = new double[N + 1];
	double *F = new double[N + 1];
	double *delta = new double[N + 1];
	//double *alpha_bl = new double[N + 1];
	double A, B, C;
	int bmax = 0;
	A = C = tau / (h_step * h_step);
	B = (1 + 2 * tau / (h_step * h_step));

	for (int i = 1; i < N; i++)	{//double *F_half = new double[N + 1];
		F[i] = (-preArray[i] - f(i * h_step, time_step));
	}
	
	delta[1] = B;
	//alpha[1] = C / B;
	beta[1] = -F[1] / B;

	for (int i = 2; i < N ; i++) {/* from 1 to N-1 */
		alpha[i-1] = C / delta[i-1];
		delta[i] = B - A * alpha[i-1];
		//alpha[i] = A / (B - C * alpha[i - i]);
		//beta[i] = (C * beta[i - 1] - F[i]) / (B - C * alpha[i - 1]);
		beta[i] = (C * beta[i - 1] - F[i]) / delta[i];
	}

	neArray[N - 1] = beta[N - 1];
	for (int i = N - 2; i > 0; i--) {
		neArray[i] = alpha[i] * neArray[i + 1] + beta[i];
	}

	delete[] delta;
	delete[] alpha;
	delete[] beta;
	delete[] F;
}

//tau/h^2 * y(i-1,j+1) - (1+2tau/h^2) * y(i,j+1) + tau/h^2*y(i+1,j+1) = -y(i,j) + tau * fi
//-y{i,j} = -y{i,j+1} + tau / h^2 *(y{i+1, j+1} - 2y{i,j+1} + y{i-1, j+1} + tau * fi
int unexplicit_scheme() {
	cout << "N=" << N << endl;
	FILE *output;
	output = fopen("unex_exact.txt", "w");
	FILE *output_1;
	output_1 = fopen("unex_solution.txt", "w");

	int k = 0;
	double max = 0;
	
	for (int i = 0; i < N + 1; i++) {
		prevd[i] = U_0(i * h_step);
		prehalfd[i] = U_0(h_step / 2 + i * h_step);
	}
	
	for (int i = 0; i < T + 1; i++) {
		nextd[0] = U_left(tau * (i + 1));
		nextd[N] = U_right(tau * (i + 1));	
		progonka(tau * i, prevd, nextd);

		nehalfd[0] = solution(h_step, tau * i + tau / 2); 
		nehalfd[N] = solution(N*h_step-h_step ,tau * i + tau / 2); //U_right(tau * i + tau / 2);
		progonka(tau * i + tau / 2, prehalfd, nehalfd);
		

		for (int j = 1; j < N; j++) {// new data become previous on next time step 
			prevd[j] = nextd[j];
			prehalfd[j] = nehalfd[j];
		}
		
		for (int j = 1; j < N; j++) {
			if (max < abs(prevd[j] - prehalfd[j]) * 2) {
				max = abs(prevd[j] - prehalfd[j]) * 2;
			}	
		}

		if (max > eps) {
			cout << max << "\nLets update the grid " << endl;
			update();
			fclose(output);
			fclose(output_1);
			return 0;
		}
		
		if (k%500 == 0) {
			for (int j = 0; j < N + 1; j++) {
				fprintf(output, "%f %f\n", j*h_step, solution(j * h_step, i * tau));
				fprintf(output_1, "%f %f\n", j*h_step, prevd[j]);
				//cout << k << " " << prevd[j] << endl;
			}
		}

		k++;
		//cout << k << endl;
				//cout << k << endl;
	}
	cout << max << endl;
	return 1;
}


int main(void) {
	int mode;
	cout << "Enter mode: 1-explicit, 2-unexplicit" << endl;
	cin >> mode;
	if (mode == 1) {
		while (explicit_scheme() == 0);
	} else {
		while (unexplicit_scheme() == 0);
	}
}
