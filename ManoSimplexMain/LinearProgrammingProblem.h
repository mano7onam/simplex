#pragma once

#include "Globals.h"
#include "Fraction.h"
#include "Matrix.h"

struct CanonicalLPP {
	int n;
	int m;
	vector<Frac> vf;
	Mat mat;
	vector<int> basis;

	CanonicalLPP(const vector<vector<Frac>> &v, vector<Frac> vf, vector<int> basis) :
		mat(v), vf(vf), basis(basis)
	{
		n = v[0].size();
		m = v.size();
		mat.subGauss(basis);
	}

	CanonicalLPP(const Mat &mat, vector<Frac> vf, vector<int> basis) :
		mat(mat), vf(vf), basis(basis)
	{
		n = mat.v[0].size();
		m = mat.v.size();
		this->mat.subGauss(basis);
	}

	CanonicalLPP(const vector<vector<Frac>> &v, vector<Frac> vf) :
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

	static CanonicalLPP getLppArtificialBasis(const CanonicalLPP &from) {
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

		CanonicalLPP res(mat, vf, basis);
		assert(res.isInCorrectBasis());
		return res;
	}

};

enum OptType {
	OPTMIN = 0, OPTMAX = 1
};
OptType getDualOptType(OptType cond) {
	return (cond == OPTMIN ? OPTMAX : OPTMIN);
}
OptType getOptTypeFromString(string sOpt) {
	if (sOpt == "min") return OPTMIN;
	else if (sOpt == "max") return OPTMAX;
	else assert(false);
}
string optTypeToString(OptType optType) {
	switch (optType) {
	case OPTMIN:
		return "min";
	case OPTMAX:
		return "max";
	default:
		assert(false);
	}
}

enum ConditionType {
	GOREQ = 0, LOREQ = 1, EQUALITY = 2
};
ConditionType getDualCondition(ConditionType cond) {
	switch (cond) {
	case GOREQ:
		return LOREQ;
	case LOREQ:
		return GOREQ;
	default:
		return EQUALITY;
	}
}
ConditionType getConditionBy(OptType opt, bool mustBePositive) {
	if (mustBePositive) {
		return (opt == OPTMIN ? GOREQ : LOREQ);
	}
	return EQUALITY;
}
ConditionType getConditionFromString(string scond) {
	if (scond == ">=") return GOREQ;
	else if (scond == "<=") return LOREQ;
	else if (scond  == "=") return EQUALITY;
	else assert(false);
}
string conditionTypeToString(ConditionType condType) {
	switch (condType) {
	case GOREQ:
		return ">=";
	case LOREQ:
		return "<=";
	case EQUALITY:
		return "=";
	default:
		assert(false);
	}
}

struct CommonLPP {
	OptType optType;
	Mat mat;
	vector<Frac> vf;
	vector<ConditionType> conditions;
	vector<bool> mustBePositive;

	CommonLPP() : mat({}) {}

	void read() {
		int n, m, k;
		cin >> n >> m >> k;

		mat.v.assign(m, vector<Frac>(n + 1, 0));
		conditions.resize(m);
		mustBePositive.assign(n, false);

		string soptType; 
		cin >> soptType; 
		optType = getOptTypeFromString(soptType);

		vf.resize(n + 1);
		for (int i = 0; i <= n; ++i) {
			vf[i].read();
		}
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mat.v[i][j].read();
			}
			string scondType;
			cin >> scondType;
			conditions[i] = getConditionFromString(scondType);

			mat.v[i][n].read();

			if (optType == OPTMIN && conditions[i] == LOREQ || optType == OPTMAX && conditions[i] == GOREQ) {
				conditions[i] = getDualCondition(conditions[i]);
				for (int j = 0; j <= n; ++j) {
					mat.v[i][j] = mat.v[i][j] * -1;
				}
			}
		}

		for (int i = 0; i < k; ++i) {
			int id;
			cin >> id;
			mustBePositive[id] = true;
		}
	}

	CommonLPP getDual() {
		CommonLPP res;
		res.mat.v.assign(mat.getN(), vector<Frac>(mat.getM() + 1, 0));
		res.optType = getDualOptType(optType);
		
		res.vf.assign(mat.getM() + 1, 0);
		for (int i = 0; i < mat.getM(); ++i) {
			res.vf[i] = mat.v[i].back();
		}
		res.vf.back() = vf.back();
		
		res.conditions.resize(mat.getN());
		for (int i = 0; i < mat.getN(); ++i) {
			for (int j = 0; j < mat.getM(); ++j) {
				res.mat.v[i][j] = mat.v[j][i];
			}
			res.mat.v[i].back() = vf[i];
			res.conditions[i] = getConditionBy(res.optType, mustBePositive[i]);
		}

		res.mustBePositive.assign(mat.getM(), false);
		for (int i = 0; i < res.mustBePositive.size(); ++i) {
			res.mustBePositive[i] = (conditions[i] != EQUALITY);
		}

		return res;
	}

	void print() {
		int n = int(mat.v[0].size()) - 1;
		int m = mat.v.size();
		int k = 0;
		for (int el : mustBePositive) k += el;
		cout << n << ' ' << m << ' ' << k << endl;
		cout << optTypeToString(optType) << endl;
		for (auto el : vf) {
			el.printA();
			cout << ' ';
		}
		cout << endl;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mat.v[i][j].printA();
				cout << ' ';
			}
			cout << conditionTypeToString(conditions[i]) << ' ';
			mat.v[i].back().printA();
			cout << endl;
		}
		for (int i = 0; i < mustBePositive.size(); ++i) {
			if (mustBePositive[i]) cout << i << ' ';
		}
	}
};
