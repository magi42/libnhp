/***************************************************************************
 *   This file is part of the NeHeP library.                               *
 *                                                                         *
 *   Copyright (C) 1997-2002 Marko Grönroos <magi@iki.fi>                  *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 *  This library is free software; you can redistribute it and/or          *
 *  modify it under the terms of the GNU Library General Public            *
 *  License as published by the Free Software Foundation; either           *
 *  version 2 of the License, or (at your option) any later version.       *
 *                                                                         *
 *  This library is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      *
 *  Library General Public License for more details.                       *
 *                                                                         *
 *  You should have received a copy of the GNU Library General Public      *
 *  License along with this library; see the file COPYING.LIB.  If         *
 *  not, write to the Free Software Foundation, Inc., 59 Temple Place      *
 *  - Suite 330, Boston, MA 02111-1307, USA.                               *
 *                                                                         *
 ***************************************************************************/

#include <ctype.h>
#include <magic/mmath.h>
#include <magic/mstream.h>
#include <magic/mtextstream.h>
#include "nhp/genetics.h"
#include "nhp/genes.h"
#include "nhp/individual.h"
#include "nhp/gaenvrnmt.h"
#include "nhp/testenv.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// -----                  ----   _   -----             o                                      //
//   |    ___   ____  |  |      / \  |       _                      _          ___    _    |  //
//   |   /   ) (     -+- | --- /   \ |---  |/ \  |   | | |/\  __  |/ \  |/|/| /   ) |/ \  -+- //
//   |   |---   \__   |  |   \ |---| |     |   |  \ /  | |   /  \ |   | | | | |---  |   |  |  //
//   |    \__  ____)   \ |___/ |   | |____ |   |   V   | |   \__/ |   | | | |  \__  |   |   \ //
////////////////////////////////////////////////////////////////////////////////////////////////

BinaryTestEAEnv::BinaryTestEAEnv (int dim, int ntargets) {
	targets.make (ntargets, dim);
	changeObjective (0);
}

void BinaryTestEAEnv::changeObjective (int target) {
	mObjective = target;

	for (int t=0; t<targets.rows; t++)
		for (int i=0; i<targets.cols; i++)
			targets.get (t, i) = mObjective;
}

void BinaryTestEAEnv::addFeaturesTo (Genome& genome) const {
	for (int i=0; i<targets.cols; i++)
		genome.add (new BinaryGene (format ("x%d", i), 1.0));
}

double dist (double x0, double y0, double x1, double y1) {
	return sqrt ((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
}

double BinaryTestEAEnv::evaluateg (const Individual& indiv) {
	double err = 0.0;
	for (int i=0; i<targets.cols; i++) {
		int x = static_cast<const BinaryGene&> (*indiv.getGene(format ("x%d", i))).getvalue();
		err += (x==targets.get (0, i))? 0.0:1.0;
		//		printf ("%g ", x);
		//		fflush (stdout);
	}
	//	printf ("%g ", err);
	//	fflush (stdout);

	return err;
}

void BinaryTestEAEnv::cycle_report (OStream& log, OStream& out) {
	;
}



//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  ----- |                -----                  ----   _   -----              //
//  |     |       ___   |    |    ___   ____  |  |      / \  |       _          //
//  |---  |  __   ___| -+-   |   /   ) (     -+- | --- /   \ |---  |/ \  |   |  //
//  |     | /  \ (   |  |    |   |---   \__   |  |   \ |---| |     |   |  \ /   //
//  |     | \__/  \__|   \   |    \__  ____)   \ |___/ |   | |____ |   |   V   O//
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

FloatTestEAEnv::FloatTestEAEnv (const StringMap& params, int d, int f) : mParams (params) {
	dim = d;
	func = f;
	mObjective = 0;
	mGeneType = ESFLOAT;
 }

void FloatTestEAEnv::addFeaturesTo (Genome& genome) const {
	for (int i=0; i<dim; i++)
		if (mGeneType==BITFLOAT)
			genome.add (new BitFloatGene (format ("x%d", i), -4.0, 4.0, 16, mParams));
		else
			genome.add (new FloatGene (format ("x%d", i), -4.0, 4.0, 1.0));
}

///////////////////////////////////////////////////////////////////////////////
//                                                     o                     //
//     |   ___   ____  |      __         _    ___   |           _    ____    //
//    -+- /   ) (     -+-    /   |   | |/ \  |   \ -+- |  __  |/ \  (        //
//     |  |---   \__   |     +-- |   | |   | |      |  | /  \ |   |  \__     //
//      \  \__  ____)   \    |    \__! |   |  \__/   \ | \__/ |   | ____)    //
//                           |                                               //
///////////////////////////////////////////////////////////////////////////////

#define pi 3.14159

#define Sum(fx)\
	double sum=0;\
	for (int i=0; i<x.size(); i++)\
		sum += fx;

#define Mul(fx)\
	double mul=1;\
	for (int i=0; i<x.size(); i++)\
		mul *= fx;

inline double sphereTF (const PackArray<double>& x) {
	Sum (sqr(x[i]));
	return sum;
}

inline double ellipsoidTF (const PackArray<double>& x) {
	Sum (sqr(double(i+1))*sqr(x[i]));
	return sum;
}

inline double negsphereTF (const PackArray<double>& x) {
	Sum (-sqr(x[i]));
	return sum;
}

inline double zerominTF (const PackArray<double>& x) {
	Sum (fabs(x[i]));
	return sum;
}

inline double TF4 (const PackArray<double>& x) {
	double k=4;
	Sum (sqr(k*x[i])-k*cos(2*pi*x[i]*k));
	return (x.size()*10 + sum)/50.0;
}

inline double TF5 (const PackArray<double>& x) {
	double k=-1000;
	Sum (-x[i]*sin(sqrt(fabs(k*x[i]))));
	return sum;
}

// Rastrigin's function
inline double TF6 (const PackArray<double>& x) {
	double k=20;
	Sum (sqr(k*4*x[i])/4000);
	Mul (cos(k*x[i]/sqrt(i+1.0)));
	return (sum-mul+1)/4.0;
}

// Generalized Rastrigin's function
inline double GenRastriginTF (const PackArray<double>& x) {
	double A=10, w=2*M_PI;
	Sum (sqr(x[i])-A*cos(w*x[i]));
	return A*x.size() + sum;
}

// Schwefel's function
inline double TF7 (const PackArray<double>& x) {
	Sum (fabs(x[i]));
	return sum;
}

// Griewangk's function
inline double TF8 (const PackArray<double>& x) {
	Sum (fabs(x[i]));
	return sum;
}


///////////////////////////////////////////////////////////////////////////////

double FloatTestEAEnv::evaluateg (const Individual& indiv) {
	PackArray<double> x (dim);
	for (int i=0; i<dim; i++)
		x[i] = dynamic_cast<const AnyFloatGene&> (*indiv.getGene(format ("x%d", i))).getvalue();
	
	return calc (x, func);
}

double FloatTestEAEnv::calc (const Vector& x0, int func) {
	ASSERT (func>=0 && func<functions);
	ASSERT (mObjective==0 || mObjective==1);

	Vector x = x0;
	if (mObjective==1)
		for (int i=0; i<x.size(); i++)
			x[i] = 1-x[i];
	
	double result=0;
	switch (func) {
	  case Sphere:		result = sphereTF		(x); break;
	  case Ellipsoid:	result = ellipsoidTF	(x); break;
	  case NegSphere:	result = negsphereTF	(x); break;
	  case ZeroMin:		result = zerominTF		(x); break;
	  case F4:			result = TF4			(x); break;
	  case F5:			result = TF5			(x); break;
	  case F6:			result = TF6			(x); break;
	  case F7:			result = TF7			(x); break;
	  case F8:			result = TF8			(x); break;
	};

	return result;
}

void FloatTestEAEnv::printMathematica2D () {
	TextOStream out;
	out.autoFlush ();

	Vector vec (2);
	for (int f=4; f<functions-2; f++) {
		out.printf ("testfunc%d := {", f);
		for (double y=-1.0; y<=1.0; y+=0.02) {
			if (y>-1.0)
				out << ",";
			out << "{";
			for (double x=-1.0; x<=1.0; x+=0.02) {
				if (x>-1.0)
					out << ",";
				vec[0] = x;
				vec[1] = y;
				double r = calc (vec, f) + gaussrnd (0.1);
				out.printf ("%f", r);
			}
			out << "}";
		}
		out.printf ("}\n"
					"ListDensityPlot[testfunc%d, Mesh -> False]\n\n",
					f);
	}
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      |   |       |     o |   | o        ----   _   -----                  //
//      |\ /|       |  |    |\ /|     _   |      / \  |       _              //
//      | V | |   | | -+- | | V | | |/ \  | --- /   \ |---  |/ \  |   |      //
//      | | | |   | |  |  | | | | | |   | |   \ |---| |     |   |  \ /       //
//      |   |  \__! |   \ | |   | | |   | |___/ |   | |____ |   |   V        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void MultiMinEAEnv::addFeaturesTo (Genome& genome) const {
	genome.add (new FloatGene ("x", 0.0, 1.0, 1.0));
	genome.add (new FloatGene ("y", 0.0, 1.0, 1.0));
}

double MultiMinEAEnv::evaluateg (const Individual& genome) {
	double	x	= static_cast<const FloatGene&> (genome["x"]).getvalue();
	double	y	= static_cast<const FloatGene&> (genome["y"]).getvalue();

	double d[4];
	d[0] = dist (x, y, 0.6, 0.75);
	d[1] = dist (x, y, 0.25, 0.25);
	d[2] = dist (x, y, 0.0, 1.0);
	d[3] = dist (x, y, 0.9, 0.4);
	
	double nearv=d[0];
	for (int i=1; i<4; i++)
		if (d[i]<nearv)
			nearv = d[i];

	return nearv;
}
