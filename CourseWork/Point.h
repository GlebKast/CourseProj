#pragma once

struct Point
{
	double x;
	double y;
	
	double alpha;
	double beta;
	double gamma;

	double fi;
	double psi;

	Point() : x(1.), y(1.), alpha(1.), beta(1.), gamma(1.), fi(1.), psi(1.) {}

	Point(double x, double y, double alpha, double beta, double gamma, double fi, double psi) :
		x(x), y(x), alpha(alpha), beta(beta), gamma(gamma), fi(fi), psi(psi) {}

	Point(Point &p) : x(p.x), y(p.y), alpha(p.alpha), beta(p.beta), gamma(p.gamma), fi(p.fi), psi(p.psi) {}

};