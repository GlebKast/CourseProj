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

	Point() : x(0.), y(0.), alpha(0.), beta(0.), gamma(0.), fi(0.), psi(0.) {}

	Point(double x, double y, double alpha, double beta, double gamma, double fi, double psi) :
		x(x), y(y), alpha(alpha), beta(beta), gamma(gamma), fi(fi), psi(psi) {}

	Point(Point &p) : x(p.x), y(p.y), alpha(p.alpha), beta(p.beta), gamma(p.gamma), fi(p.fi), psi(p.psi) {}

};