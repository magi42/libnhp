/***************************************************************************
 *   This file is part of the NeHeP library.                               *
 *                                                                         *
 *   Copyright (C) 1997-2005 Marko Grönroos <magi@iki.fi>                  *
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

/*
  Strategian Tien perusta kouluni näkökulmasta selitetään Maan
  kirjassa. On vaikeaa ymmärtää todellinen Tie vain miekkailun
  avulla. Tiedä pienimmät asiat ja suurimmat asiat, matalimmat asiat
  ja syvimmät asiat. Kuin se olisi suora tie kirjoitettuna maahan,
  ensimmäistä kirjaa kutsutaan Maan kirjaksi.
  */

#ifndef __EAENVRNMNT_H__
#define __EAENVRNMNT_H__

#include <magic/mobject.h>
#include <magic/mexception.h>

using namespace MagiC;

/** Reference type. EXPERIMENTAL.
 **/
/*
template <class TYPE>
class Ref : public Object {
  protected:
	TYPE*	obj;
  public:
				Ref						(TYPE* o) {obj = o;}
				Ref						(Ref<TYPE>& o) {obj = o.obj; o.obj=NULL;}
				~Ref					() {delete obj;}
				operator =				(TYPE* n) {delete obj; obj=n;}
	//				operator const TYPE&	() const {return *obj;}
				operator TYPE&			() {return *obj;}
	TYPE*		operator ->				() {return obj;}
	TYPE*		getptr					() {return obj;}
	//	TYPE*		operator &				() {return obj;}
	//	const TYPE*	operator &				() const {return obj;}
};
*/

// Externals
class Genome;
class Individual;

/** The abstract base class for environments ("objective functions")
 *  where the fitness of Individuals is measured.
 *
 *  The inheritors of this class are typically used as follows (assume
 *  that the parameter map has been defined and parameters to it
 *  inserted):
 *
 *	MyEnvironment myEnvir (paramMap, someOtherParameters);
 *
 *  SimplePopulation pop (envir, paramMap);
 *
 *	pop.evolve (100);
 *  
 *  Note: the class is actually not abstract because the RTTI system
 *  doesn't allow that.
 **/
class EAEnvironment : public Object {
	decl_dynamic (EAEnvironment);
  public:

					EAEnvironment	();

	/** Sets the number of evaluations to be done for each individual
	 *  to determine it's fitness. This is useful only if the fitness
	 *  measurement is noisy.
	 *
	 *  This variable should be set before starting the evolution, of
	 *  course.
	 **/
	void			setevals		(int n) {mNEvals = n;}

	/** Returns the number of evaluations to be done for each
	 *  individual to determine it's fitness.
	 **/
	int				evals			() const {return mNEvals;}

	/** Returns the total number of evaluations performed in all
	 *  generations so far. This includes all multiple evaluations of
	 *  all individuals.
	 **/
	int				total_evals		() const {return mTotEvals;}

	/** Orders to add gaussian noise to the evaluated fitness. This is
	 *  useful for testing the robustness of the EA for noisy
	 *  evaluation.
	 **/
	void			addnoise		(double stddev) {mNoise = stddev;}

	double			evaluate		(const Individual& indiv);

	/** Sets evolution log directory for cycle reports and
	 *  miscellaneous log files.
	 **/
	void			logDir			(const String& dir) {mLogDir = dir;}

	/** Prints out a generation report to the given brief log stream
	 *  and more verbose output stream.
	 **/
	void			cycleReport		(OStream& log, OStream& out) {mCycles++; cycle_report(log,out);}
	
	// Virtual methods
	
	/** Adds problem-specific features to a genome.
	 **/
	virtual void	addFeaturesTo	(Genome& genome) const {;}

	/** Initializes a generation. Not necessary in all models.
	 **/
	virtual void	init_cycle		() {mBestFitness = 99999999.0; mpBest = NULL;}

	/** Implementation for @ref Object. */
	virtual void		check			() const;
	
	// Lowest non-subjective fitness; unlike the fitness returned by
	// the evaluate-method, this does NOT include the artificial noise
	// factor.
	double			bestfitn;

	// ThreadLock&		getLock() const {return mLock;}
	
  protected:
	/** Evaluates the fitness of the given individual in the
	 *  environment. MUST OVERLOAD!
	 *
	 *  @exception must_overload
	 **/
	virtual double	evaluateg		(const Individual& ind) {MUST_OVERLOAD; return 0.0;}

	/** Prints some statistics or something at the end of the evaluation cycle.
	 **/
	virtual void	cycle_report	(OStream& log, OStream& out) {;}

	int				mNEvals; //< Number of evaluations to be done for each individual, default=1.
	int				mTotEvals; //< Total number of evaluations done so far.
	int				mCycles; //< Basicly a generation counter. It would be nicer if this was not stored here (redundacy, see).

	/** Amount of artificial noise to be added to evaluated
	 *  fitness. This is useful for testing the effects of noise on
	 *  evolution. Requires deterministic test functions to be useful.
	 **/
	double			mNoise;
	Individual*		mpBest; //< The best individual in this cycle.

	/** (Subjective) fitness of the best individual. "Subjective"
	 *  means here that the artificial noise is included in this value.
	 **/
	double			mBestFitness;
	String			mLogDir; //< Directory for gathering evolution logs.

	// mutable ThreadLock	mLock;
};

#endif
