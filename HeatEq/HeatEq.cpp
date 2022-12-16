#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>

using namespace std;

//Зададим константы
const double pi = 3.14159265;
double eps = 0.01;

double sigma = 1.0;
double betta = sigma + 1;
double k0 = 0.001;
double q0 = 0.0005;
double tf = 100.0;

double T = 50.0; //Граница по времени
double tau = 0.1; //Шаг по времени

double L_T = 2 * pi * sqrt(k0 / q0) * sqrt((sigma + 1) / sigma * sigma);
double l = L_T; //Границы где точное решение >=0
double step = 0.1; //шаг по пространству

int N1 = (2 * L_T / step); //количество шагов по пространству
int N2 = (T / tau); //Количество шагов по времени

//Минимальный и максимальный шаг по времени для контроля точности
double tau_min = 0.001; 
double tau_max = 1;

double gamma = tau / (step * step);

double* coord, * y_null, * y_curr, * y_next, * alpha, * beta, * delta, * exact, * A, * B, * C, * F;

/*источник тепла*/
double f(double u)
{
	return q0 * powl(u, betta);
}

/* коэф теплопроводности нелинейный*/
double k(double u)
{
	//return u;
	return k0 * powl(u, sigma);
}

double u0(double x) /* U(x,0) = U_0(x)*/
{
	return 2 * sin(x); // В данном случае не используется, так как начальное распределение считаем по страшной формуле из Самарского Тихонова 
}

// граничные условия
double u1(double t) /* U(0,t) = U_1(t)*/
{
	return 0.0;
}

double u2(double t) /* U(l,t) = U_2(t)*/
{
	return 0.0;
}

double max2(double* y1, double* y2)
{
	double max = 0;

	for (int i = 0; i < N1 + 1; i++)
		if (abs(y1[i] - y2[i]) > max)
			max = abs(y1[i] - y2[i]);
	
	return max;
}

double max1(double* y1)
{
	double max = 0;

	for (int i = 0; i < N1 + 1; i++)
		if (abs(y1[i]) > max)
			max = abs(y1[i]);

	return max;
}

//Для рассчета погрешности
double error(double* y1, double* y2)
{
	double error = 0;

	for (int i = 0; i < N1 + 1; i++)
		error = error + (y1[i] - y2[i]) * (y1[i] - y2[i]);

	error = sqrt(error);

	return error;
}

//Прогонка
void one_step(double* y, double tau, double time, double* new_y)
{
	double accuracy = 0;
	int count = 0;
	gamma = tau / (step * step);

	double* y_s1, * y_s;

	y_s1 = (double*)calloc(N1 + 1, sizeof(double)); //s1 операция
	y_s = (double*)calloc(N1 + 1, sizeof(double)); //s операция

	for (int i = 0; i < N1 + 1; i++)
		y_s1[i] = y[i];

	for (int i = 0; i < N1 + 1; i++)
		y_curr[i] = y[i];

	do
	{
		count++;

		for (int i = 0; i < N1 + 1; i++)
			y_s[i] = y_s1[i];

		A[0] = 0;
		A[1] = 0;
		for (int i = 2; i < N1; i++)
			A[i] = (k0 * gamma / 2) * (pow(y_s[i - 1], sigma) + pow(y_s[i], sigma));

		C[0] = 0;
		for (int i = 1; i < N1; i++)
			C[i] = -1 - (k0 * gamma / 2) * (pow(y_s[i + 1], sigma) + 2 * pow(y_s[i], sigma) + pow(y_s[i - 1], sigma));

		B[0] = 0;
		for (int i = 1; i < N1 - 1; i++)
			B[i] = (k0 * gamma / 2) * (pow(y_s[i + 1], sigma) + pow(y_s[i], sigma));

		for (int i = 0; i < N1; i++)
			F[i] = -y_curr[i] - tau * q0 * pow(y_s[i], betta);

		alpha[0] = 0; alpha[1] = 0; alpha[2] = -B[1] / C[1];
		beta[0] = 0; beta[1] = 0;
		beta[1] = F[1] / C[1];


		for (int i = 3; i < N1; i++)
		{
			alpha[i] = -B[i - 1] / (A[i - 1] * alpha[i - 1] + C[i - 1]);
			beta[i] = (F[i - 1] - A[i - 1] * beta[i - 1]) / (A[i - 1] * alpha[i - 1] + C[i - 1]);
		}

		//Не забыть про граничные условия
		y_s1[0] = u1(time);
		y_s1[N1] = u2(time);

		y_s1[N1 - 1] = (F[N1 - 1] - A[N1 - 1] * beta[N1 - 1]) / (C[N1 - 1] + A[N1 - 1] * alpha[N1 - 1]);
		for (int i = N1 - 2; i > 0; i--)
			y_s1[i] = alpha[i + 1] * y_s1[i + 1] + beta[i + 1];

		//Проконстрлируем погрешность
		accuracy = max2(y_s1, y_s) / max1(y_s1);
		if ((accuracy <= eps) || (count == 2))
			break;

	} while (true);

	//Скопируем массив
	for (int i = 0; i < N1 + 1; i++)
		new_y[i] = y_s1[i];

}

//Неявная схема методом прогонки
void implicit_scheme(double* y, double tau)
{
	double current_time = 0.0;
	FILE* Local_out, * Local_time, * Local_final, * Local_exact;
	fopen_s(&Local_out, "Implicit_Localization_out.txt", "w");
	fopen_s(&Local_time, "Implicit_Localization_tau.txt", "w");
	fopen_s(&Local_final, "Implicit_Localization_final.txt", "w");
	fopen_s(&Local_exact, "Localization_exact.txt", "w");

	double* y_tau1 = (double*)calloc(N1 + 1, sizeof(double)); //одиночны шаг
	double* y_tau2 = (double*)calloc(N1 + 1, sizeof(double)); //половинный шаг

	int cnt = 0;
	do
	{
		cnt++;

		if (current_time >= T)
			break;

		one_step(y, tau, current_time + tau, y_tau1);
		one_step(y, tau / 2, current_time + tau / 2, y_tau2);
		one_step(y, tau / 2, current_time + tau / 2, y_tau2);

		//Контроль точности
		if (error(y_tau1, y_tau2) < eps)
		{
			do
			{
				tau = tau * 2;
				one_step(y, tau, current_time + tau, y_tau1); /* tau/2 -> tau */
				one_step(y, tau / 2, current_time + tau / 2, y_tau2); /* tau/2 -> tau*/
				one_step(y_tau2, tau / 2, current_time + tau / 2, y_tau2);
			} while ((error(y_tau1, y_tau2) < eps) & (tau_min <= tau) & (tau <= tau_max));

			tau = tau / 2;
		}
		else
		{
			do
			{
				tau = tau / 2;
				one_step(y, tau, current_time + tau, y_tau1); /* tau *2 - > tau */
				one_step(y, tau / 2, current_time + tau / 2, y_tau2); /* tau -> tau/2 */
				one_step(y_tau2, tau / 2, current_time + tau / 2, y_tau2);
			} while ((error(y_tau1, y_tau2) > eps) & (tau_min <= tau) & (tau <= tau_max));

			tau = tau * 2;
		}

		//cout << "tau=" << tau << endl;
		cout << "current_time=" << current_time << endl;
		current_time = current_time + tau;
		one_step(y, tau, current_time, y);

		if (cnt % 10 == 0)
			for (int i = 0; i < N1 + 1; i++)
				fprintf(Local_out, "%f %f\n", coord[i], y[i]);

		fprintf(Local_time, "%f \n", tau);

		for (int i = 0; i < N1 + 1; i++)
		{
			double koef = 2 * (sigma + 1) / (sigma * (sigma + 2));
			double cosin = cos((coord[i] * pi) / L_T) * cos((coord[i] * pi) / L_T);
			exact[i] = powl(q0 * (tf - T), (-1 / sigma)) * powl(koef * cosin, 1 / sigma);
			if (abs(coord[i]) > l / 2)
				exact[i] = 0;
			fprintf(Local_exact, "%f %f\n", coord[i], exact[i]);
		}

	} while (current_time < T);

	for (int i = 0; i < N1 + 1; i++)
		fprintf(Local_final, "%f %f\n", coord[i], y[i]);

	std::fclose(Local_out);
	std::fclose(Local_time);
}

int main()
{
	FILE* Local_start;
	fopen_s(&Local_start, "Localization_start.txt", "w");

	//Сетка по пространству
	coord = (double*)calloc(N1 + 1, sizeof(double));

	//Начальное распределение
	y_null = (double*)calloc(N1 + 1, sizeof(double));

	//Данные на i-ом слое
	y_curr = (double*)calloc(N1 + 1, sizeof(double));
	
	//Данные на i+1 слое
	y_next = (double*)calloc(N1 + 1, sizeof(double)); /* data on j+1 layer */

	//Точное решение
	exact = (double*)calloc(N1 + 1, sizeof(double)); /* array for exact */

	//Для прогонки
	alpha = (double*)calloc(N1, sizeof(double));
	beta = (double*)calloc(N1, sizeof(double));
	delta = (double*)calloc(N1, sizeof(double));
	A = (double*)calloc(N1, sizeof(double));
	B = (double*)calloc(N1, sizeof(double));
	C = (double*)calloc(N1, sizeof(double));
	F = (double*)calloc(N1 + 1, sizeof(double));

	//cout << "Coord step " << N1 << endl;
	//cout << "Time step " << N2 << endl;
	//cout << "L= " << l << endl;

	//Зададим сетку
	for (int i = 0; i < N1 + 1; i++)
		coord[i] = l - i * step;

	double cosin;
	double koef;
	//Начальное распределение
	for (int i = 0; i < N1 + 1; i++)
	{
		koef = 2 * (sigma + 1) / (sigma * (sigma + 2));
		cosin = cos((coord[i] * pi) / L_T) * cos((coord[i] * pi) / L_T);
		y_null[i] = powl(q0 * tf, (-1 / sigma)) * powl(koef * cosin, 1 / sigma);
		if (abs(coord[i]) > l / 2)
			y_null[i] = 0;
		fprintf(Local_start, "%f %f\n", coord[i], y_null[i]);
	}

	implicit_scheme(y_null, tau);

	return 0;
}