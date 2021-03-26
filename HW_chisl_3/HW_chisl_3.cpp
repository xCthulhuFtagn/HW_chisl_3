#include <vector>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <cmath>
#include <exception>
#include <algorithm>
#include <windows.h>
#include <conio.h>
#undef min
#undef max
using namespace std;
/*
vector<double>& vector<double>::operator=(const vector<double>& right) {
	this->assign(right.front(), right.back());
}
void operator-=(vector<double>& left, const vector<double>& right) {
	if (left.size() != right.size()) throw exception("Differently sized vectors!");
	else {
		unsigned size = left.size();
		for (auto i = 0; i < size; ++i) {
			left[i] -= right[i];
		}
	}
}*/

vector<double> operator*(const vector<double>& left, const double right) {
	vector<double> ans;
	for (auto el : left) {
		ans.push_back(el * right);
	}
	return ans;
}

double Lagrange(valarray<double>* data, double param) {
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

double Newton(valarray<double>* data, double param) {
	unsigned size = data[0].size(), i;
	vector<vector<double>> f(size - 1);
	double ans, tmp;
	for (i = 0; i < size - 1; ++i) {
		f[0].push_back((data[1][i + 1] - data[1][i]) / (data[0][i + 1] - data[0][i]));
	}
	for (unsigned i = 1; i < size - 1; ++i) {
		for (unsigned j = 0; j < size - i - 1; ++j) {
			f[i].push_back((f[i - 1][j + 1] - f[i - 1][j]) / (data[0][j + i + 1] - data[0][j]));
		}
	}
	ans = data[1][0];
	for (unsigned i = 0; i < size - 1; ++i) {
		tmp = f[i].front();
		for (unsigned j = 0; j <= i; ++j) {
			tmp *= (param - data[0][j]);
		}
		ans += tmp;
	}
	return ans;
}

valarray<double> Gaus(vector <valarray<double>> a) {
	unsigned n = a.size();
	double coef;
	valarray<double> ans(n);
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
		ans[i]=(a[i][n] / a[i][i]);
	}
	return ans;
}

double Spline(valarray<double>* input, double param) {
	if (param <= input[0][0]) throw invalid_argument("Can't find Splines for such argument!");
	vector<double> h;
	double ans;
	unsigned size = input[0].size(), i;
	vector<valarray<double>> output(size);
	valarray<double> m;
	for (unsigned i = 0; i < size - 1; ++i) {
		h.push_back(input[0][i + 1] - input[0][i]);
	}
	for (i = 0; i < size - 2; ++i) {
		output[i].resize(size + 1, 0);
		output[i][i + 1] = (h[i] + h[i + 1]) / 3;
		output[i][i + 2] = h[i + 1] / 6;
		output[i][size] = (input[1][i + 2] - input[1][i + 1]) / h[i + 1] - (input[1][i + 1] - input[1][i]) / h[i];
	}
	output[i].resize(size + 1, 0);
	output[i][0] = 1;
	output[i + 1].resize(size + 1, 0);
	output[i + 1][size - 1] = 1;
	m = Gaus(output);
	for (i = 0; i < size && input[0][i] < param; ++i);
	ans = ((m[i] * pow((param - input[0][i - 1]), 3) + m[i - 1] *
		pow(input[0][i] - param, 3)) / (h[i - 1] * 6) + (input[1][i] - m[i] * pow(h[i - 1], 2) / 6) * ((param - input[0][i - 1]) / h[i - 1]) + (input[1][i - 1] - m[i - 1] * pow(h[i - 1], 2) / 6) * (input[0][i] - param) / h[i - 1]);
	return ans;
}

double SquareAproximation(valarray<double>* data, double x) {
	unsigned size = data[0].size();
	vector<valarray<double>> suem(3);
	valarray<double> abc;
	vector<double> koef(7, 0);
	for (unsigned i = 0; i < 7; ++i) {
		if (i < 4) {
			for (unsigned j = 0; j < size; ++j) {
				koef[i] += pow(data[0][j], i + 1);//иксы
			}
		}
		else {
			for (unsigned j = 0; j < size; ++j) {
				koef[i] += data[1][j] * pow(data[0][j], i - 4);//игреки
			}
		}
	}
	for (unsigned i = 0; i < 3; ++i) {
		suem[i].resize(4);
		for (unsigned j = 0; j < 3; ++j) {
			if (i != 2 || j != 2) {
				suem[i][j] = koef[3 - j - i];
			}
			else {
				suem[i][j] = size;
			}
		}
		suem[i][3] = koef[6 - i];
	}
	abc = Gaus(suem);
	return abc[0] * pow(x, 2) + abc[1] * x + abc[2];
}

double LinearApproximation(valarray<double>* data, double x)
{
	unsigned size = data[0].size();
	vector<double> koef(4, 0);
	valarray<double> ans;
	for (int i = 0; i < 4; i++) {
		koef[0] += data[0][i];
		koef[1] += pow(data[0][i], 2);
		koef[2] += data[1][i];
		koef[3] += data[0][i] * data[1][i];
	}
	vector<valarray<double>> syst(2);
	for (auto& el : syst) el.resize(3);
	syst[0][0] = koef[1];
	syst[0][1] = koef[0];
	syst[0][2] = koef[3];
	syst[1][0] = koef[1];
	syst[1][1] = size;
	syst[1][2] = koef[2];
	ans = Gaus(syst);
	return ans[0] * x + ans[1];
	//тут надо еще посмотреть
	double Discrepancy = 0;
	for (int i = 1; i <= size; i++) {
		Discrepancy += pow(ans[0] * i + ans[1] - syst[i][2], 2);
	}
}

void Draw(valarray<double>* data, double(*f)(valarray<double>*, double)) {
	Sleep(1000);
	printf("\033[2J");
	double x_max, x_min, diap_x;
	double y_max, y_min, diap_y;
	double y, x;
	y_min = data[1].min(); x_min = data[0].min();
	y_max = data[1].max(); x_max = data[0].max();
	diap_x = (x_max - x_min) / 10;
	diap_y = (y_max - y_min) / 10;                             //диапазон делений для оY 
	y_max += diap_y; y_min -= diap_y;
	HDC hDC = GetDC(GetConsoleWindow());                       //настройка консоли для рисования 
	HPEN Pen = CreatePen(PS_SOLID, 1, RGB(255, 255, 255));     //ручка для разметки 
	HPEN Pen1 = CreatePen(PS_SOLID, 2, RGB(0, 255, 0));        //ручка для графика 
	HPEN Pen2 = CreatePen(PS_SOLID, 8, RGB(255, 0, 0));        //ручка для точек
	if (true) {
		HWND hwnd = GetConsoleWindow();
		RECT rect;
		//через прямоугольник rect описывается консолька 
		GetWindowRect(hwnd, &rect);
		int width = rect.right - rect.left;
		int height = rect.bottom - rect.top;
		SelectObject(hDC, Pen);
		//оси координат 
		MoveToEx(hDC, 0, height / 2, NULL);
		LineTo(hDC, width, height / 2);
		MoveToEx(hDC, width / 2, height, NULL);
		LineTo(hDC, width / 2, 0);
		//стрелка oY 
		MoveToEx(hDC, width / 2, 0, NULL);
		LineTo(hDC, width / 2 + 10, 10);
		MoveToEx(hDC, width / 2 - 10, 10, NULL);
		LineTo(hDC, width / 2, 0);
		//стрелка оХ 
		MoveToEx(hDC, width - 45, height / 2 - 10, NULL);
		LineTo(hDC, width - 35, height / 2);
		MoveToEx(hDC, width - 35, height / 2, NULL);
		LineTo(hDC, width - 45, height / 2 + 10);
		//деления oX
		for (x = x_min; x <= x_max; x += diap_x) {
			MoveToEx(hDC, width * (x - x_min) / (x_max - x_min), height / 2 - height / 250 - 1, NULL);
			LineTo(hDC, width * (x - x_min) / (x_max - x_min), height / 2 + height / 250 + 1);
		}
		//деления oY
		for (y = y_min; y < y_max + diap_y; y += diap_y) {
			MoveToEx(hDC, width / 2 - width / 250 - 1, height * (y_max - y) / (y_max - y_min), NULL);
			LineTo(hDC, width / 2 + width / 250 + 1, height * (y_max - y) / (y_max - y_min));
		}
		//вывод графика и точек, по которому он строится
		SelectObject(hDC, Pen1);
		for (x = x_min; x < x_max; x += 1e-3) {
			MoveToEx(hDC, width * (x - x_min) / (x_max - x_min), height * (y_max - f(data, x)) / (y_max - y_min), NULL);
			LineTo(hDC, width * (x + 1e-3 - x_min) / (x_max - x_min), height * (y_max - f(data, x + 1e-3)) / (y_max - y_min));
		}
		SelectObject(hDC, Pen2);
		for (unsigned i = 0; i < data[0].size(); ++i) {
			MoveToEx(hDC, width * (data[0][i] - x_min) / (x_max - x_min), height * (y_max - data[1][i]) / (y_max - y_min), NULL);
			LineTo(hDC, width * (data[0][i] - x_min) / (x_max - x_min), height * (y_max - data[1][i]) / (y_max - y_min));
		}
	}
	_getch();
}


int main()
{
	unsigned n, size;
	double LagPol = 0, x, y, param, tmp;
	valarray<double> data[2], splines;
	try {
		cout << "Enter number of input: ";
		cin >> n;
		if (n == 0 || !cin) throw invalid_argument("Wrong number");
		data[0].resize(n);
		data[1].resize(n);
		cout << "Enter data:" << endl;
		for (unsigned i = 0; i < n; ++i) {
			cin >> x >> y;
			if (!cin) throw invalid_argument("Wrong data");
			data[0][i]=(x);
			data[1][i]=(y);
		}
		cout << "Choose what you want to see:\nInterpolation by\n1)Newton\n2)Lagrange\n3)Splines\nApproximations :\n4)Quadratic\n5)Linear\n";
		cin >> n;
		if (!cin || n <= 0 || n > 4) throw invalid_argument("Wrong number entered while choosing what to see");
		else {
			switch (n) {
			case 1:
				Draw(data, Newton);
			case 2:
				Draw(data, Lagrange);
			case 3:
				Draw(data, SquareAproximation);
			case 4:
				Draw(data, LinearApproximation);
			}
		}
	}
	catch (invalid_argument& bad) {
		cout << bad.what();
		return -1;
	}
	return 0;
}

