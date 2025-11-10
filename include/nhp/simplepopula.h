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

#ifndef __SIMPLEPOPULA_H__
#define __SIMPLEPOPULA_H__

#include "population.h"

#include <magic/mthread.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  ---- o            |       ----                  |           o            //
// (              --  |  ___  |   )       --        |  ___   |           _   //
//  ---  | |/|/| |  ) | /   ) |---   __  |  ) |   | |  ___| -+- |  __  |/ \  //
//     ) | | | | |--  | |---  |     /  \ |--  |   | | (   |  |  | /  \ |   | //
// ___/  | | | | |    |  \__  |     \__/ |     \__! |  \__|   \ | \__/ |   | //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** A linear set implementation of @ref Population.
 **/
class SimplePopulation : public Population {
  public:

	/** Standard constructor.
	 **/
								SimplePopulation (EAEnvironment& envr, const StringMap& params);

								~SimplePopulation ();

	/** Evolves the population.
	 *
	 *  @param generations Number of generations the population should be evolved.
	 *  @param savefile File where the evolution log should be saved.
	 *
	 *  @param trg_fitn The evolution is terminated if the attained
	 *  fitness goes smaller than this value. Special value -1 tells
	 *  not to use this rule.
	 **/
	double						evolve			(int generations,
												 const char* savefile=NULL,
												 double trg_fitn=-1);
	
	/** Returns an individual in the population by it's index number.
	 **/
	const Individual&			operator[]		(int i) const {return (*mpPopulation)[i];}

	/** Returns an individual in the population by it's index
	 *  number. As above, except non-conest.
	 **/
	Individual&					operator[]		(int i) {return (*mpPopulation)[i];}
	
	/** Returns the number of individuals in the population.
	 **/
	int							size			() const {return mpPopulation->size();}

	/** Dumps the population to given output in a formatted manner.
	 **/
	void						print			(TextOStream& out) const;

	/** Writes an one-line generation report to the given logging stream.
	 **/
	void						report			(TextOStream& log) const;

	/** Returns the current strategy.
	 **/
	const EAStrategy&			getstrategy		() const {return *mpStrategy;}

	/** Returns the number of times the population has been evaluated
	 *  (i.e. the generation)
	**/
	int							getAge			() const {return mAge;}

	/** Returns a non-const reference to the global mutation rate.
	 **/
	MutationRate&				mutRate			() {return mGlobalMutationRate;}

	/** Returns the selection parameters.
	 **/
	const SelectionPrms&		selParams		() const {return mSelectionParams;}
	/** Returns a non-const reference to the selection parameters.
	 **/
	SelectionPrms&				selParams		() {return mSelectionParams;}

	/** Resets the stored fitness averages of multiply measured
	 *  individuals. Useful when changing the objective function (old
	 *  fitness values would be invalid).
	 **/
	void						resetFitnesses	();

	/** Implementation for @ref Object. */
	virtual void				check			() const;
	
	/** Minimum allowed similarity between genomes.
	 **/
	double						minsimilarity;

	const Array<Individual>&	getPopArray	() const {return *mpPopulation;}
  private:

	/** Implementation for @ref Population. Evaluates the population
	 *  in the given environment.
	 **/
	void					evaluate		(EAEnvironment& envr, TextOStream& out);
	void					evaluate		(int i, EAEnvironment& environment, TextOStream& out);
		
 	/** Implementation for @ref Population. Add population-dependent
	 *  features to a genome. I suppose there might be some use for
	 *  this. Maybe.
	 **/
	void					addFeaturesTo	(Genome& genome) const;

  private:
	Array<Individual>*		mpPopulation;		/**> The current set of individuals in the population. */
	EAStrategy*				mpStrategy;			/**> The evolutionary strategy, the evolutinary algorithms. */
	FitnessStats			mFitnessStats;		/**> Some statistics about the current population. */
	int						mElites;			/**> Number of elites, individuals who should survive intact to the next generation. */
	bool					mUseGlobalElites;	/**> Mode flag indicating whether or not elites should be used. */
	SelectionPrms			mSelectionParams;	/**> Selection parameters. */
	bool					mUseGlobalMu;		/**> Global portion of the number of potential parents (mu). The semantics are dependent on the evolutionary strategy used. */
	bool					mUseGlobalQ;		/**> Global q parameter for tournament selection. */
	bool					mUseGlobalEtaPlus;	/**> Global etaPlus parameter for linear ranking selection. */
	ThreadLock				mThreadLock;        /**> For locking data. */
	
	friend class EAStrategy;
	friend class Selector;
	friend class EvaluationWorker;
};

#endif
