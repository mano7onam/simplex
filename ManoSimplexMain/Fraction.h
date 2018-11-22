#pragma once

#include "Globals.h"

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

	void printA() {
		cout << a;
	}

	void read2() {
		cin >> a >> b;
	}
};
