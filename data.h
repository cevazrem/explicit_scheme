#pragma once

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

void init()
{
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
}

void update()
{
	if (flag == 1)
	{
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
	else
	{
		flag = 2;
	}
}
