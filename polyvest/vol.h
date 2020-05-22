/***********************************************************************
*  This code is part of PolyVest.
*
*  Copyright (C) 2013, 2016 Cunjing Ge, Institute of Software, Chinese 
*  Academy of Sciences, Beijing, China. All rights reserved. 
*  E-mail: <gecj@ios.ac.cn>.
*
*  PolyVest is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
* 
*  PolyVest is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with PolyVest. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include <iostream>
#include <algorithm>
#include <cmath>
#include "armadillo"
#include "time.h"
#include "memory.h"

#ifndef POLYVOL_H
#define POLYVOL_H

using namespace std;
using namespace arma;

namespace vol{

class Polyvest_p{
public:
	Polyvest_p(int rows, int cols);
	~Polyvest_p();
	void 	matA(double val, int i, int j){ A(i, j) = val; }
	double 	matA(int i, int j){ return A(i, j); }
	void 	vecb(double val, int i){	b(i) = val; }
	double 	vecb(int i){ return b(i); }
	void 	genInitE(double &R2, vec &Ori);
	void 	Preprocess_hacked();
	void 	Preprocess();
	double 	EstimateVol(int coef);
	double 	Volume() const { return vol; }
	mat A;
	vec b;

	bool 	msg_off;
	bool 	check_planes_off;
private:
	double 	walk(int k);
	void	checkHPs();

	double 	randd(double u){ return rand() * u / RAND_MAX; }
	int 	randi(int u){ return rand() % u; }

//polytope denoted by: Ax<=b, A is an (m x n) matrix.
	int m, n;

//approximating volume variables
	vec x;
	double vol, determinant;
	int l;
	double *r2;

	//reciprocal of beta, beta-cut
	double beta_r;

	vec *B;
	mat *Ai;
};

inline Polyvest_p::Polyvest_p(int rows, int cols) :
    msg_off(true),
	check_planes_off(false),
	m(rows),
	n(cols),
	A(rows, cols),
	b(rows),
	x(cols),
	vol(0),
	determinant(0)
{
	srand((int)time(0));

	beta_r = 2 * n; //2 * n;

	l = (int)(n * log((double)beta_r) / log((double)2)) + 2;
	r2 = new double[l];
	for (int i = 0; i < l; i++) 
		r2[i] = pow((double)2, (double)(2 * i) / n);

	B = new vec[n];
	Ai = new mat[n];
}

inline Polyvest_p::~Polyvest_p(){
	delete []r2;
	delete []B;
	delete []Ai;
}
	
}

#endif
