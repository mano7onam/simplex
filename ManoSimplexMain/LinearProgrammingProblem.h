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

struct CommonLPP {
	Mat mat;
	vector<Frac> vf;

};
