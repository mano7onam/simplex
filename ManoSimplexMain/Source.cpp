#pragma warning (disable : 4996)

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cmath>
#include <sstream>

using namespace std;

#define mp make_pair
typedef long long ll;
typedef long double ld;

ll gcd(ll a, ll b) {
	return b == 0 ? a : gcd(b, a % b);
}

int llLen(ll num) {
	stringstream ss;
	ss << num;
	string s;
	ss >> s;
	return s.size();
}

struct Frac {
	ll a, b;

	int getLen() {
		return llLen(a) + 1 + llLen(b) + 2;
	}

	void normalize() {
		ll g = gcd(abs(a), abs(b));
		a /= g;
		b /= g;
	}

	void notNegB() {
		if (b < 0) {
			b = -b;
			a = -a;
		}
	}

	Frac(ll a = 0, ll b = 1) : a(a), b(b) {
		notNegB();
		assert(this->b > 0);
		normalize();
	}

	void read() {
		cin >> a;
		b = 1;
	}

	Frac rev() const { 
		assert(a != 0);
		return Frac(b, a); 
	}

	bool operator == (const Frac &f) const { return a * f.b == b * f.a; }
	bool operator != (const Frac &f) const { return !(operator==(f)); }
	bool operator < (const Frac &f) const { return a * f.b < b * f.a; }
	bool operator >= (const Frac &f) const { return !(operator<(f)); }
	bool operator <= (const Frac &f) const { return (operator<(f)) || (operator==(f)); }

	Frac operator + (const Frac& f) const { return Frac(a * f.b + f.a * b, b * f.b); }
	Frac operator - (const Frac& f) const { return Frac(a * f.b - f.a * b, b * f.b); }
	Frac operator * (const Frac& f) const { return Frac(a * f.a, b * f.b); }
	Frac operator / (const Frac& f) const { 
		return operator*(f.rev());
	}

	void print(int len = 10) {
		if (0 == a) {
			cout << "(0)";
			for (int cur = 3; cur < len; ++cur) {
				cout << ' ';
			}
			return;
		}
		cout << "(" << a << "/" << b << ")";
		for (int cur = getLen(); cur < len; ++cur) {
			cout << ' ';
		}
	}
};

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
			for (int i = row + 1; i < v.size(); ++i) {
				if (v[i][col] == 0) continue;
				v[i] = vsub(v[i], vmul(v[row], v[i][col]));
			}
			print();
		}
		row++;
		for (int ic = int(cols.size()) - 1; ic >= 0; ic--) {
			row--;
			int col = cols[ic];
			int rowNext = row + 1;
			for (int icNext = ic + 1; icNext < cols.size(); ++icNext) {
				int colNext = cols[icNext];
				if (v[row][colNext] == 0) continue;
				v[row] = vsub(v[row], vmul(v[rowNext], v[row][colNext]));
				rowNext++;
			}
			print();
		}
	}
};

struct SimplexTable {
	vector<vector<Frac>> v;
	vector<int> basis;

	void print_basis() {
		cout << "Basis: ";
		for (auto el : basis) cout << el << ' ';
		cout << endl;
	}

	void print() {
		cout << "----------- Simplex table ----------------" << endl;
		for (int i = 0; i < v.size(); ++i) {
			for (int j = 0; j < v[i].size(); ++j) {
				v[i][j].print();
				cout << ' ';
			}
			cout << endl;
		}
		print_basis();
		cout << "------------------------------------------" << endl;
		cout << endl;
	}

	SimplexTable(int n, int m, vector<Frac> vf, const Mat &mat, const vector<int> &cols) {
		n++; m++;
		v.assign(m, vector<Frac>(n, 0));
		basis = cols;

		for (int i = 0; i < mat.v.size(); ++i) {
			for (int j = 0; j < mat.v[i].size(); ++j) {
				int nj = (j + 1) % n;
				v[i + 1][nj] = mat.v[i][j];
			}
		}

		for (int row = 0; row < cols.size(); ++row) {
			int col = cols[row];
			auto koef = vf[col];
			vf[col] = 0;
			vf = vsub(vf, vmul(mat.v[row], koef));
			vf[col] = 0;
		}

		for (int j = 0; j < vf.size(); ++j) {
			int nj = (j + 1) % n;
			v[0][nj] = vf[j];
		}

		print();
	}

	void doSimplexMethod() {
		for (int i = 1; i < v.size(); ++i) {
			assert(v[i][0] >= 0);
		}

		bool haveAns = false;
		while (true) {
			bool good = true;
			for (int s = 1; s < v[0].size(); ++s) {
				if (v[0][s] < 0) {
					good = false;
					break;
				}
			}
			if (good) {
				haveAns = true;
				break;
			}

			int bCol = -1, bRow = -1;
			Frac bRes;
			for (int s = 1; s < v[0].size(); ++s) {
				if (v[0][s] >= 0) continue;
				for (int r = 1; r < v.size(); ++r) {
					if (v[r][s] <= 0) continue;
					Frac curRes = v[r][0] / v[r][s];
					if (-1 == bCol || curRes < bRes) {
						bRes = curRes;
						bCol = s;
						bRow = r;
					}
				}
			}
			if (-1 == bCol) break;

			basis[bRow - 1] = bCol - 1;

			for (int r = 0; r < v.size(); ++r) {
				if (r == bRow) {
					v[r] = vmul(v[r], v[bRow][bCol].rev());
				}
				else {
					auto koef = v[r][bCol] / v[bRow][bCol];
					v[r] = vsub(v[r], vmul(v[bRow], koef));
				}
			}

			print();
		}

		if (haveAns) {
			auto res = v[0][0] * (-1);
			cout << "Have answer:" << endl;
			cout << "w = ";
			res.print();
			cout << endl;
			cout << "Basis: ";
			for (auto el : basis) cout << el << ' ';
			cout << endl;
			cout << "Solution: ( ";
			vector<Frac> ans(v[0].size() - 1);
			for (int i = 0; i < basis.size(); ++i) {
				ans[basis[i]] = v[i + 1][0];
			}
			for (int i = 0; i < ans.size(); ++i) {
				ans[i].print();
				if (i < int(ans.size()) - 1) cout << ", ";
			}
			cout << " )";
			cout << endl;
		}
		else {
			cout << "No answer." << endl;
		}
	}
};

void solve() {
	int n, m;
	cin >> n >> m;
	vector<Frac> vf(n + 1);
	for (int i = 0; i < vf.size(); ++i) {
		vf[i].read();
	}
	vector<vector<Frac>> v(m, vector<Frac>(n + 1, 0));
	for (int i = 0; i < v.size(); ++i) {
		for (int j = 0; j < v[i].size(); ++j) {
			v[i][j].read();
		}
	}
	vector<int> cols(m);
	for (int i = 0; i < m; ++i) {
		cin >> cols[i];
	}
	Mat mat(v);
	mat.subGauss(cols);

	SimplexTable st(n, m, vf, mat, cols);
	st.doSimplexMethod();
}

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	solve();
	// just for commit

	return 0;
}