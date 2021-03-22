#include <vector>
#include <iostream>
#include <iomanip>
#include <valarray>
using namespace std;

vector<double>& vector<double>::operator=(const vector<double>& right){
	this->assign(right.front(), right.back());
}

double Lagrange(double param, vector<double>* data) {
	double tmp, LagPol = 0;
	unsigned n = data[0].size();
	for (unsigned i = 0; i < n; ++i) {
		tmp = data[1][i];
		for (unsigned j = 0; j < n; ++j) {
			if (i == j) continue;
			tmp *= (param - data[0][j]) / (data[0][i] - data[0][j]);
		}
		LagPol += tmp;
	}
	return LagPol;
}

double Newton(double param, vector<double>* data) {
	unsigned size = data[0].size(), i;
	vector<vector<double>> f(size-1);
	double ans, tmp;
	for (i = 0; i < size - 1; ++i) {
		f[0].push_back((data[1][i + 1] - data[1][i]) / (data[0][i + 1] - data[0][i]));
	}
	for (unsigned i = 1; i < size-1; ++i) {
		for (unsigned j = 0; j < size - i -1; ++j) {
			f[i].push_back((f[i - 1][j + 1] - f[i - 1][j]) / (data[0][j + i + 1] - data[0][j])); 
		}
	}
	ans = data[1].front();
	for (unsigned i = 0; i < size-1; ++i) {
		tmp = f[i].front();
		for (unsigned j = 0; j <= i; ++j) {
			tmp *= (param - data[0][j]);
		}
		ans += tmp;
	}
	return ans;
}

vector<double> Gaus(vector <valarray<double>> a) {
	unsigned n = a.size();
	double coef;
	vector<double> ans;
	for (unsigned i = 0; i < n; ++i) {
		unsigned j;
		for (j = i; j < n && a[j][i] == 0; ++j);
		if (j == n) {
			cout << "System has no solution or their quantity is infinity" << endl;
			exit(0);
		}
		else if (i != j)
			swap(a[j], a[i]);
		for (j = 0; j < i; ++j) {
			coef = a[j][i] / a[i][i];
			a[j] -= a[i] * coef;
		}
		for (j = i + 1; j < n; ++j) {
			coef = a[j][i] / a[i][i];
			a[j] -= a[i] * coef;
		}
	}
	for (unsigned i = 0; i < n; ++i) {
		ans.push_back(a[i][n] / a[i][i]);
	}
	return ans;
}

double Spline(double param, vector<double> *input) {
	if (param <= input[0][0]) throw invalid_argument("can't find Splines for such argument!");
	vector<double> h, m;
	double ans;
	unsigned size = input[0].size(), i;
	vector<valarray<double>> output(size);
	for (unsigned i = 0; i < size-1; ++i) {
		h.push_back(input[0][i + 1] - input[0][i]);
	}
	for (i = 0; i < size-2; ++i) {
		output[i].resize(size + 1, 0);
		output[i][i + 1] = (h[i] + h[i + 1]) / 3;
		output[i][i + 2] = h[i + 1] / 6;
		output[i][size] = (input[1][i + 2] - input[1][i + 1]) / h[i + 1] - (input[1][i + 1] - input[1][i]) / h[i];
	}
	output[i].resize(size + 1, 0);
	output[i][0] = 1;
	output[i+1].resize(size + 1, 0);
	output[i+1][size-1] = 1;
	m=Gaus(output);
	for (i = 0; i < size && input[0][i] < param; ++i);
	ans=((m[i] * pow((param - input[0][i-2]), 3) + m[i-1] * pow(input[0][i-1] - param, 3)) / (h[i-1] * 6) + (input[1][i-1] - m[i] * pow(h[i-1], 2) / 6) * ((param - input[0][i-2]) / h[i-1]) + (input[1][i-2] - m[i-1] * pow(h[i-1], 2) / 6) * (input[0][i-1] - param) / h[i-1]);
	return ans;
}

int main()
{
	unsigned n, size;
	double LagPol = 0, x, y, param, tmp;
	vector<double> data[2], splines;
	try {
		cout << "Enter number of input: ";
		cin >> n;
		if (n == 0 || !cin) throw invalid_argument("wrong number");
		cout << "Enter data:" << endl;
		for (unsigned i = 0; i < n; ++i) {
			cin >> x >> y;
			if (!cin) throw invalid_argument("wrong data");
			data[0].push_back(x);
			data[1].push_back(y);
		}
	}
	catch (invalid_argument bad) {
		cout << bad.what();
		return -1;
	}
	size = data[0].size();
	cout << "Lagrange:" << endl;
	cout << "x\t|\ty\t|\tans" << endl;
	for (unsigned i = 0; i < size; ++i) {
		cout << data[0][i] << "\t|\t" << data[1][i] << "\t|\t" << Lagrange(data[0][i], data) << endl;
	}
	cout << "Newton:" << endl;
	cout << "x\t|\ty\t|\tans" << endl;
	for (unsigned i = 0; i < size; ++i) {
		cout << data[0][i] << "\t|\t" << data[1][i] << "\t|\t" << Newton(data[0][i], data) << endl;
	}
	cout << "Bias in 0.1875 is: " << Lagrange(0.1875, data) - Newton(0.1875, data) << endl;
	try {
		cout << Spline(0.1875, data);
	}
	catch (invalid_argument& bad) {
		cout << bad.what();
	}
	return 0;
}