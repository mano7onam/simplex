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

enum ConditionType {
	GOREQ = 0, LOREQ = 1, EQUALITY = 2
};
ConditionType getDualCondition(ConditionType cond) {
	if (cond == GOREQ) {
		return LOREQ;
	}
	else if (cond == LOREQ) {
		return GOREQ;
	}
	return EQUALITY;
}
ConditionType getConditionBy(OptType opt, bool mustBePositive) {
	if (mustBePositive) {
		return (opt == OPTMIN ? GOREQ : LOREQ);
	}
	return EQUALITY;
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

		mat.v.assign(m, vector<Frac>(n, 0));
		conditions.resize(m);
		mustBePositive.assign(n, false);

		string soptType; 
		cin >> soptType; 
		if (soptType == "min") optType = OPTMIN;
		else if (soptType == "max") optType = OPTMAX;
		else assert(false);

		vf.resize(n + 1);
		for (int i = 0; i < n; ++i) {
			vf[i].read();
		}
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mat.v[i][j].read();
			}
			string scondType;
			cin >> scondType;
			if (scondType == ">=") conditions[i] = GOREQ;
			else if (scondType == "<=") conditions[i] = LOREQ;
			else if (scondType == "=") conditions[i] = EQUALITY;
			else assert(false);

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
};
