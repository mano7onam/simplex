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

const ll INF = 1e9 + 7;

ll gcd(ll a, ll b) {
	return b == 0 ? a : gcd(b, a % b);
}

struct Frac {
	ll a, b;

	int llLen(ll num) {
		stringstream ss;
		ss << num;
		string s;
		ss >> s;
		return s.size();
	}

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

struct LPP {
	int n;
	int m;
	vector<Frac> vf;
	Mat mat;
	vector<int> basis;

	LPP (const vector<vector<Frac>> &v, vector<Frac> vf, vector<int> basis) :
		mat(v), vf(vf), basis(basis)
	{
		n = v[0].size();
		m = v.size();
		mat.subGauss(basis);
	}

	LPP(const Mat &mat, vector<Frac> vf, vector<int> basis) :
		mat(mat), vf(vf), basis(basis)
	{
		n = mat.v[0].size();
		m = mat.v.size();
		this->mat.subGauss(basis);
	}

	LPP (const vector<vector<Frac>> &v, vector<Frac> vf) :
		mat(v), vf(vf)
	{
		n = v[0].size();
		m = v.size();
		basis = vector<int>(m);
		for (int i = 0; i < basis.size(); ++i) basis[i] = i;
		mat.subGauss(basis);
	}

	bool isInCorrectBasis() const {
		for (int i = 0; i < m; ++i) {
			if (mat.v[i].back() < 0) return false;
		}
		return true;
	}

	void changeBasis(const vector<int> &newBasis) {
		basis = newBasis;
		mat.subGauss(basis);
	}

	static LPP getLppArtificialBasis(const LPP &from) {
		assert(!from.isInCorrectBasis());
		int mini = -1;
		Frac minv(INF, 1);
		for (int i = 0; i < from.m; ++i) {
			if (-1 == mini || from.mat.v[i].back() < minv) {
				minv = from.mat.v[i].back();
				mini = i;
			}
		}
		assert(minv < 0);

		int n = from.n + 1;
		vector<Frac> vf(n, 0);
		vf[n - 2] = 1;
		
		auto mat = from.mat;
		int m = from.m;
		auto basis = from.basis;
		for (int i = 0; i < m; ++i) {
			mat.v[i].push_back(-1);
			int sz = int(mat.v[i].size());
			swap(mat.v[i][sz - 1], mat.v[i][sz - 2]);
			if (mini == i) {
				basis[i] = sz - 2;
			}
		}

		LPP res(mat, vf, basis);
		assert(res.isInCorrectBasis());
		return res;
	}

};

struct SimplexTable {
	vector<vector<Frac>> v;
	vector<int> basis;
	bool haveAns;

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

	SimplexTable(const LPP &lpp) {
		haveAns = false;
		v.assign(lpp.m + 1, vector<Frac>(lpp.n, 0));
		int height = int(v.size());
		int width = int(v[0].size());
		basis = lpp.basis;
		auto &mat = lpp.mat;
		auto vf = lpp.vf;

		for (int i = 0; i < mat.v.size(); ++i) {
			for (int j = 0; j < mat.v[i].size(); ++j) {
				int nj = (j + 1) % width;
				v[i + 1][nj] = mat.v[i][j];
			}
		}

		for (int row = 0; row < basis.size(); ++row) {
			int col = basis[row];
			auto koef = vf[col];
			vf[col] = 0;
			vf = vsub(vf, vmul(mat.v[row], koef));
			vf[col] = 0;
		}

		for (int j = 0; j < vf.size(); ++j) {
			int nj = (j + 1) % width;
			v[0][nj] = vf[j];
		}

		print();
	}

	vector<int> getBasisIndices() {
		return basis;
	}

	vector<Frac> getBasisVector() {
		vector<Frac> ans(v[0].size() - 1);
		for (int i = 0; i < basis.size(); ++i) {
			ans[basis[i]] = v[i + 1][0];
		}
		return ans;
	}

	Frac getAns() {
		return v[0][0] * (-1);
	}

	Mat getMatrix() {
		int height = int(v.size());
		int width = int(v[0].size());
		vector<vector<Frac>> res(height - 1, vector<Frac>(width, 0));
		for (int i = 1; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				int ni = i - 1;
				int nj = (j + width - 1) % width;
				res[ni][nj] = v[i][j];
			}
		}
		return res;
	}

	void printAns() {
		if (haveAns) {
			auto res = getAns();
			cout << "Have answer:" << endl;
			cout << "w = ";
			res.print();
			cout << endl;
			cout << "Basis: ";
			for (auto el : basis) cout << el << ' ';
			cout << endl;
			cout << "Solution: ( ";
			auto ans = getBasisVector();
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

	void doSimplexMethod() {
		for (int i = 1; i < v.size(); ++i) {
			assert(v[i][0] >= 0);
		}

		haveAns = false;
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
	vector<int> basis(m);
	for (int i = 0; i < m; ++i) {
		cin >> basis[i];
	}

	LPP lpp(v, vf, basis);
	if (lpp.isInCorrectBasis()) {
		SimplexTable st(lpp);
		st.doSimplexMethod();
		st.printAns();
		cout << "It was received with 1 phase\n";
	}
	else {
		auto sublpp = LPP::getLppArtificialBasis(lpp);
		SimplexTable subst(sublpp);
		subst.doSimplexMethod();
		subst.printAns();
		if (subst.getAns() != 0) {
			cout << "No solution. Can not to find basic valid solution\n";
			return;
		}
		cout << "\n\n-------------------------------------------------------\n";
		cout << "-------------------------------------------------------\n";
		cout << "-------------------------------------------------------\n\n";

		auto mat = subst.getMatrix();
		int sz = int(mat.v[0].size());
		basis = subst.getBasisIndices();
		for (int row = 0; row < basis.size(); ++row) {
			int b = basis[row];
			if (b == sz - 2) {
				for (int j = 0; j < sz; ++j) {
					if (mat.v[row][j] != 0) {
						basis[row] = j;
						break;
					}
				}
				assert(basis[row] < sz - 2);
				break;
			}
		}

		lpp.changeBasis(basis);

		SimplexTable st(lpp);
		st.doSimplexMethod();
		st.printAns();
		cout << "It was received with 2 phases\n";
	}
}

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	solve();

	return 0;
}