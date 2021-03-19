#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

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

int main()
{
	unsigned n, size;
	double LagPol = 0, x, y, param, tmp;
	vector<double> data[2];
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
	cout << "x\t|\ty\t|\tans" << endl;
	for (unsigned i = 0; i < size; ++i) {
		cout << data[0][i] << "\t|\t" << data[1][i] << "\t|\t" << Lagrange(data[0][i], data) << endl;
	}
	return 0;
}