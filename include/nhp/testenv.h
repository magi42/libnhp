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

#include "nhp/gaenvrnmt.h"
#include <magic/mtable.h>
#include <magic/mmath.h>
#include <magic/mmap.h>

////////////////////////////////////////////////////////////////////////////////////////////////
// -----                  ----   _   -----             o                                      //
//   |    ___   ____  |  |      / \  |       _                      _          ___    _    |  //
//   |   /   ) (     -+- | --- /   \ |---  |/ \  |   | | |/\  __  |/ \  |/|/| /   ) |/ \  -+- //
//   |   |---   \__   |  |   \ |---| |     |   |  \ /  | |   /  \ |   | | | | |---  |   |  |  //
//   |    \__  ____)   \ |___/ |   | |____ |   |   V   | |   \__/ |   | | | |  \__  |   |   \ //
////////////////////////////////////////////////////////////////////////////////////////////////

/** Environment for making EA experiments with a binary-valued genomes.
 **/
class BinaryTestEAEnv : public EAEnvironment {
	// Should be bool-table but gcc doesn't like them
	PackTable<int>	targets;
	int				mObjective;
  public:

	/** Creates a binary test environment.
	 *
	 *  @param d Dimension, i.e. number of binary genes in the genome.
	 *
	 *  @param nObjectives Number of alternative minima (objectives)
	 *  in the search problem.
	**/
					BinaryTestEAEnv	(int d, int nObjectives=1);

	/** Changes objective to another target (usually 0 or 1).
	 **/
	void			changeObjective				(int target);
	
	// Implementations
	
	virtual void	addFeaturesTo				(Genome& genome) const;
	virtual void	init_cycle					() {;}
	virtual double	evaluateg					(const Individual& genome);
	virtual void	cycle_report				(OStream& log, OStream& out);
};



//////////////////////////////////////////////////////////////////////////////////
//  ----- |                -----                  ----   _   -----              //
//  |     |       ___   |    |    ___   ____  |  |      / \  |       _          //
//  |---  |  __   ___| -+-   |   /   ) (     -+- | --- /   \ |---  |/ \  |   |  //
//  |     | /  \ (   |  |    |   |---   \__   |  |   \ |---| |     |   |  \ /   //
//  |     | \__/  \__|   \   |    \__  ____)   \ |___/ |   | |____ |   |   V   O//
//////////////////////////////////////////////////////////////////////////////////

/** Environment for testing @ref FloatGene and @ref BitFloatGene genes.
 **/
class FloatTestEAEnv : public EAEnvironment {
  public:

	/** Standard constructor.
	 *
	 *  @param params Dynamic parameters in a @ref String @ref
	 *  Map. Important only for defining ["BitFloatGene.grayCoding"].
	 *
	 *  @param dim Dimension of the search space, i.e. number of
	 *  floating-point genes in the genome. Genes will gave value
	 *  range [-4,4]. BitFloatGenes will have 16 bits.
	 *
	 *  @param funct_b ?
	**/
					FloatTestEAEnv		(const StringMap& params,
											 int dim=2, int funct_b=-1);

	/** Makes an extensive search of all search spaces (functions) in
	 *  two dimensions. The output is printed as Mathematica matrices.
	 **/
	void			printMathematica2D	();

	/** Calculates the value of test function with given input vector.
	 *
	 *  @param v Parameter vector.
	 *  @param f ID of test function.
	 **/
	double	calc						(const Vector& v, int f);

	/** Changes the objective.
	 **/
	void			changeObjective		(int o) {mObjective=o;}

	/** Sets the gene type: name ESFLOAT is @ref FloatGene and name
	 *  BITFLOAT is @ref BitFloatGene (with 16 bits).
	 **/
	void			setGeneType			(int vt) {mGeneType=vt;}
	
	// Implementations
	
	virtual void	addFeaturesTo			(Genome& genome) const;
	virtual void	init_cycle				() {;}
	virtual double	evaluateg				(const Individual& genome);
	virtual void	cycle_report			(OStream& log, OStream& out) {;}

	enum testfunctions {Sphere=0, Ellipsoid, NegSphere, ZeroMin,
						F4, F5, F6, F7, F8, /* Number of functions: */ functions};
	enum genetypes {ESFLOAT=0, BITFLOAT};

  protected:
	/** Dimension of search space. */
	int dim;

	/** Function. */
	int	func;

	/** Function-dependent parameter. */
	double par_a;

	/** Function-dependent parameter. */
	double par_b;

	int		mObjective;
	int		mGeneType;

	const StringMap& mParams;
	
};

/** @ref EAEnvironment for testing two-dimensional multiple-minima test function. INCOMPLETE!
 **/
class MultiMinEAEnv : public EAEnvironment {
  public:

	/** Creates a test environment. minima is the number of
	 *  minima. minimadist is the desired distance between
	 *  these minima.
	 *
	 *  @param minima Number of minima
	 *  @param minimadist Distance between the minima.
	 **/
					MultiMinEAEnv	(int minima, int minimadist) {;}
	virtual void	addFeaturesTo	(Genome& genome) const;
	virtual void	init_cycle		() {;}
	virtual double	evaluateg		(const Individual& genome);
	virtual void	cycle_report	(OStream& log, OStream& out) {;}
};


