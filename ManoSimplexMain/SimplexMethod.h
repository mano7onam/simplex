#pragma once

#include "Globals.h"
#include "Matrix.h"
#include "LinearProgrammingProblem.h"

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

	SimplexTable(const CanonicalLPP &lpp) {
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

void readAndSolveLinearProgrammingProblemCanonical(bool readInitialBasis) {
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

	CanonicalLPP lpp(v, vf);

	if (readInitialBasis) {
		vector<int> basis(m);
		for (int i = 0; i < m; ++i) {
			cin >> basis[i];
		}
		lpp.basis = basis;
	}
	
	if (lpp.isInCorrectBasis()) {
		SimplexTable st(lpp);
		st.doSimplexMethod();
		st.printAns();
		cout << "It was received with 1 phase\n";
	}
	else {
		auto sublpp = CanonicalLPP::getLppArtificialBasis(lpp);
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
		auto basis = subst.getBasisIndices();
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
