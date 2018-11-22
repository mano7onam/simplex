#pragma warning (disable : 4996)

#include "Globals.h"
#include "SimplexMethod.h"
#include "LinearProgrammingProblem.h"

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	//readAndSolveLinearProgrammingProblemCanonical(false);
	CommonLPP lpp;
	lpp.read();
	CommonLPP dualLpp = lpp.getDual();
	dualLpp.print();

	return 0;
}