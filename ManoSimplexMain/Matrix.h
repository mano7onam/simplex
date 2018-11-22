#pragma once

#include "Globals.h"
#include "Fraction.h"

vector<Frac> vmul(const vector<Frac> &v, Frac koef) {
	vector<Frac> res;
	for (auto el : v) {
		res.push_back(el * koef);
	}
	return res;
}

vector<Frac> vadd(const vector<Frac> &a, const vector<Frac> &b) {
	vector<Frac> res;
	for (int i = 0; i < a.size(); ++i) {
		res.push_back(a[i] + b[i]);
	}
	return res;
}

vector<Frac> vsub(const vector<Frac> &a, const vector<Frac> &b) {
	vector<Frac> res;
	for (int i = 0; i < a.size(); ++i) {
		res.push_back(a[i] - b[i]);
	}
	return res;
}

struct Mat {
	vector<vector<Frac>> v;

	Mat(const vector<vector<Frac>> &v) : v(v) {}

	int getN() {
		return v[0].size() - 1;
	}

	int getM() {
		return v.size();
	}

	void print() {
		cout << "---------- Matrix ---------------" << endl;
		for (int i = 0; i < v.size(); ++i) {
			for (int j = 0; j < v[i].size(); ++j) {
				v[i][j].print();
				cout << ' ';
			}
			cout << endl;
		}
		cout << "---------------------------------" << endl;
		cout << endl;
	}

	void subGauss(vector<int> cols) {
		int row = -1;
		for (int col : cols) {
			row++;
			for (int i = row; i < v.size(); ++i) {
				if (v[i][col] != 0) {
					swap(v[row], v[i]);
					break;
				}
			}
			if (v[row][col] == 0) continue;
			v[row] = vmul(v[row], v[row][col].rev());
			bool changed = false;
			for (int i = row + 1; i < v.size(); ++i) {
				if (v[i][col] == 0) continue;
				changed = true;
				v[i] = vsub(v[i], vmul(v[row], v[i][col]));
			}
			if (changed) print();
		}
		row++;
		for (int ic = int(cols.size()) - 1; ic >= 0; ic--) {
			row--;
			int col = cols[ic];
			int rowNext = row + 1;
			bool changed = false;
			for (int icNext = ic + 1; icNext < cols.size(); ++icNext) {
				int colNext = cols[icNext];
				if (v[row][colNext] == 0) continue;
				changed = true;
				v[row] = vsub(v[row], vmul(v[rowNext], v[row][colNext]));
				rowNext++;
			}
			if (changed) print();
		}
	}
};
