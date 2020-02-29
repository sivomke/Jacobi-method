#include <iostream>
#include<iomanip>
using namespace std;

const int n = 10;
const double  eps = 0.00000001;
double sign(double x) {
	if (x > 0.) return 1.;
	else return -1.;
}

double aux_sum(double(*S)[n], double(*D)[n], int i, int j) {
	double sum = 0;
	for (int l = 0; l < i; ++l) {
		sum += S[l][i] * S[l][j] * D[l][l];
	}
	return sum;
}
void Func(double(*A)[n], double(*S)[n], double(*D)[n]) {
	for (int i = 0; i < n; ++i) {
		D[i][i] = sign(A[i][i] - aux_sum(S, D, i, i));
		S[i][i] = sqrt(abs(A[i][i] - aux_sum(S, D, i, i)));
		for (int j = 0; j < n; ++j) {
			if (i > j) {
				S[i][j] = 0;
			}
			else S[i][j] = (A[i][j] - aux_sum(S, D, i, j)) / (S[i][i] * D[i][i]);
		}
	}
}

void traspose(double(*S)[n], double(*S_tr)[n]) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			S_tr[i][j] = S[j][i];
		}
	}
}

void aux_mult(double(*S)[n], double(*D)[n], double(*C)[n]) {
	for (int j = 0; j < n; ++j) {
		for (int i = j; i < n; ++i) {
			C[i][j] = S[j][i] * D[j][j];
		}
	}
}

void up_down(double(*C)[n], double* b, double* y) {
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			sum += C[i][j] * y[j];
		}
		y[i] = (b[i] - sum) / C[i][i];
		sum = 0;
	}
}

void down_up(double(*S)[n], double* y, double* x) {
	double sum = 0;
	for (int i = n - 1; i >= 0; --i) {
		for (int j = i + 1; j < n; ++j) {
			sum += S[i][j] * x[j];
		}
		x[i] = (y[i] - sum) / S[i][i];
		sum = 0;
	}
}

double det(double(*S)[n], double(*D)[n]) {
	double prod = 1;
	for (int i = 0; i < n; ++i) {
		prod *= D[i][i] * S[i][i] * S[i][i];
	}
	return prod;
}

double** inverse_matrix(double(*C)[n], double(*S)[n]) {
	double** inverse = new double* [n];
	for (int i = 0; i < n; ++i) {
		double e[n] = { 0 };
		double y[n] = { 0 };
		double* x = new double[n];
		e[i] = 1;
		up_down(C, e, y);
		down_up(S, y, x);
		inverse[i] = x;
	}
	return inverse;
}

void print_matrix(double(*A)[n]) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << setw(15) << right << fixed << A[i][j];
		}
		std::cout << endl;
	}
	std::cout << endl;
}

void print_vector(double* x) {
	for (int i = 0; i < n; ++i) {
		std::cout << setw(15) << right << fixed << x[i];
	}
}

double vector_norm(double* x) {
	double max = abs(x[0]);
	for (int i = 1; i < n; ++i) {
		if (abs(x[i]) > max) max = abs(x[i]);
	}
	return max;
}

double matrix_norm(double(*A)[n]) {
	double sum = 0, max = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			sum += abs(A[i][j]);
		}
		if (sum > max) {
			max = sum;
		}
		sum = 0;
	}
	return max;
}

double cond(double(*A)[n], double** inverse_A) {
	double sum = 0, max = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			sum += abs(inverse_A[i][j]);
		}
		if (sum > max) {
			max = sum;
		}
		sum = 0;
	}
	double inverse_A_norm = max;
	return matrix_norm(A) * inverse_A_norm;
}

double** identity_matrix(double(*A)[n], double** inverse_A) {
	double** E = new double* [n];
	double prod_sum = 0;
	for (int k = 0; k < n; ++k) {
		double* x = new double[n];
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				prod_sum += A[i][j] * inverse_A[j][k];
			}
			x[i] = prod_sum;
			prod_sum = 0;
		}
		E[k] = x;
	}
	return E;
}

double error_norm(double(*A)[n], double* x, double* b) {
	double b_calc[n] = { 0 };
	cout << "b calculated: " << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			b_calc[i] += A[i][j] * x[j];
		}
		cout << setw(15) << right << fixed << b_calc[i];
	}
	cout << endl;
	double eps[n] = { 0 };
	cout << "Error vector: " << endl;
	for (int i = 0; i < n; ++i) {
		eps[i] = b[i] - b_calc[i];
		cout << setw(15) << right << scientific << eps[i];
	}
	cout << endl;
	return vector_norm(eps);
}



/*Jacoby*/

void Jacoby(double(*A)[n], double* x, double* b) {
	double x_cur[n] = { 0 };
	double error[n] = { 0 };
	double error_ = 0;
	int iter = 0;
	do {
		for (int i = 0; i < n; ++i) {
			x_cur[i] = b[i];
			for (int j = 0; j < n; ++j) {
				if (j != i) {
					x_cur[i] -= A[i][j] * x[j];
				}
			}
			x_cur[i] /= A[i][i];
			error[i] = x_cur[i] - x[i];
		}
		error_ = vector_norm(error);
		for (int i = 0; i < n; ++i) {
			x[i] = x_cur[i];
		}
		iter++;
	} while (error_ > eps);
	cout << setw(15) << right << "Iterations: " << iter << endl;
}



int main()
{
	int i = 0, j = 0;
	std::cout.precision(8);
	double A[n][n] = { 0 };
	double b[n] = { 0 };
	double y[n] = { 0 };
	double x[n] = { 0 };
	double S[n][n] = { 0 };
	double S_tr[n][n] = { 0 };
	double C[n][n] = { 0 };
	double D[n][n] = { 0 };


	for (i = 1; i < n + 1; i++) {
		A[i - 1][i - 1] = 1. / i + i * n * n;
		for (j = i + 1; j < n + 1; j++) {
			A[i - 1][j - 1] = pow(-1, j) * j;
			A[j - 1][i - 1] = A[i - 1][j - 1];
		}
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			S_tr[i][j] = S[j][i];
		}
	}

	cout << "Matrix A:" << endl;
	print_matrix(A);
	cout << endl;
	for (int i = 1; i < n + 1; ++i) {
		b[i - 1] = 1. / i;
	}
	cout << "Vector b:" << endl;
	print_vector(b);
	cout << endl;
	Func(A, S, D);
	cout << "Matrix D:" << endl;
	print_matrix(D);
	cout << endl;
	cout << "Matrix S:" << endl;
	print_matrix(S);

	std::cout << endl;

	traspose(S, S_tr);
	cout << "Matrix S transposed:" << endl;
	print_matrix(S_tr);
	cout << endl;
	std::cout << endl;

	aux_mult(S, D, C);
	cout << "Matrix C:" << endl;
	print_matrix(C);

	std::cout << endl;

	up_down(C, b, y);
	cout << "Vector y:" << endl;
	print_vector(y);

	std::cout << endl;

	down_up(S, y, x);
	cout << "SOLUTION:" << endl;
	print_vector(x);

	std::cout << endl;

	std::cout << "Det: " << det(S, D) << endl;
	std::cout << endl;
	cout << "Inverse matrix A:" << endl;
	double** inverse_A = inverse_matrix(C, S);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(15) << right << fixed << inverse_A[i][j];
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << setw(15) << right << fixed << "Matrix A norm: " << matrix_norm(A) << endl;

	cout << endl;


	cout << setw(15) << right << fixed << "Matrix A cond: " << cond(A, inverse_A) << endl;

	cout << "Product of A and inverse A:" << endl;
	double** identity = identity_matrix(A, inverse_A);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(15) << right << fixed << identity[i][j];
		}
		cout << endl;
	}
	cout << "Error norm:" << setw(15) << right << scientific << error_norm(A, x, b) << endl;

	double guess[n] = { 0 };
	cout << "Jacoby method:" << endl;
	Jacoby(A, guess, b);
	cout << endl;
	cout << endl;
	cout << "SOLUTION:" << endl;
	print_vector(guess);
	cout << endl;
	cout << "Error norm:" << setw(15) << right << scientific << error_norm(A, x, b) << endl;
	delete[] identity;
	delete[] inverse_A;
	system("pause");
}
