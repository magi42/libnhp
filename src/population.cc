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

#include <magic/mmath.h>
#include <magic/mstream.h>

#include "nhp/population.h"
#include "nhp/simplepopula.h"
#include "nhp/gaenvrnmt.h"
#include "nhp/selection.h"
#include "nhp/mutrecord.h"

// For mutrecord.h
bool MutabilityRecord::record=false;		// Should we record or not
int MutabilityRecord::smFloatSamples=0;
int MutabilityRecord::smFloatVarSamples=0;
int MutabilityRecord::smBoolSamples=0;
double MutabilityRecord::smFloatVarMin=0;
double MutabilityRecord::smFloatVarSum=0;
double MutabilityRecord::smFloatVarMax=0;
double MutabilityRecord::smFloatMin=0;
double MutabilityRecord::smFloatSum=0;
double MutabilityRecord::smFloatMax=0;
double MutabilityRecord::smBoolSum=0;
double MutabilityRecord::smBoolMin=0;
double MutabilityRecord::smBoolMax=0;


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           -----   _    ----                                               //
//           |      / \  (      |       ___   |   ___                        //
//           |---  /   \  ---  -+- |/\  ___| -+- /   )  ___  \   |           //
//           |     |---|     )  |  |   (   |  |  |---  (   \  \  |           //
//           |____ |   | ___/    \ |    \__|   \  \__   ---/   \_/           //
//                                                      __/   \_/            //
///////////////////////////////////////////////////////////////////////////////

EAStrategy::EAStrategy (SimplePopulation& pop) : mrPopula (pop) {
	allow_same_parents = false;
	// selmethod = NULL;

	//
	// As a speed optimization we use two populations so that we can
	// use copy() operations for creating a descendant instead of
	// cloning them each time
	// 
	// By cloning the population each time
	// 14.20u 1.80s 0:16.29 98.2%
	// By copying the population to a double-buffered corpse population
	// 9.21u 1.87s 0:11.08 100.0%

	// Create the corpse vector
	mpNextGen = new Array<Individual> ();
	mpNextGen->make (mrPopula.size());

	// Create corpses for the future descendants.
	// As you may notice, there is _always_ kept some empty room at
	// the beginning of the mpNextGen to place possible elites there
	for (int i=mrPopula.mElites; i<mrPopula.size(); i++)
		mpNextGen->put (new Individual (mrPopula[0]), i);

}

void EAStrategy::evolve (EAEnvironment& envr, TextOStream& out, TextOStream& log)
{
	FUNCTION_BEGIN;
	
	envr.init_cycle ();

	// Evaluate
	mrPopula.evaluate (envr, out);
	
	// Print out cycle reports
	// out << "Reporting...\n";
	mrPopula.report (log);
	envr.cycleReport (log, out);

	SelectionSituation situation (mrPopula);
	// Order by fitness. Selection methods can use this order if they wish
	//RefArray<Individual> pop_order (*mrPopula.mpPopulation);
	//pop_order.quicksort ();
	
	// Create a selection matrix
	SelectionMatrix selmat (situation);

	// Create the next generation according to the selection matrix
	recombine (situation, selmat);

	// Re-evaluate Tarzan a little...
	if (mrPopula.mElites>0) {
		// out << "Re-evaluating Tarzan...\n";
		Individual& tarzan = const_cast<Individual&> (situation.getOrdered(0));
		tarzan.evaluate (envr, true);
		tarzan.addking ();
	}

	FUNCTION_END;
}

/*******************************************************************************
 *
 ******************************************************************************/
void EAStrategy::recombine (
	const SelectionSituation& situation,
	const SelectionMatrix&    selmat)
{
	////////////////////////////////////////////////////////////////////////////
	// Clone the old population. This is done in two parts because
	// we don't want to replicate the elites for no reason

	// Copy the elite references. Warning! Elites are now owned by two
	// objects for a while.
	for (int i=0; i<mrPopula.mElites; i++)
		mpNextGen->put (situation.getOrdered(i), i);

	////////////////////////////////////////////////////////////////////////////
	// Recombine
	
	const Individual *parent_a, *parent_b;
	for (int i=mrPopula.mElites; i<mrPopula.size(); i++) {
		// Select two parents
		int parent_a_ind, parent_b_ind;
		selmat.selectRandomPair (parent_a_ind, parent_b_ind);
		parent_a = &situation.getOrdered (parent_a_ind);
		parent_b = &situation.getOrdered (parent_b_ind);
		
		// Recombine them as the descendant
		(*mpNextGen)[i].recombine (*parent_a, *parent_b);

		// Mutate the descendant a little
		(*mpNextGen)[i].pointMutate (mrPopula.mutRate());

		// Incarnate the descendant
		(*mpNextGen)[i].incarnate (true);
	}

	////////////////////////////////////////////////////////////////////////////
	// Finally, set the new population as current

	// Swap the next generation as current
	Array<Individual>* tmp = mrPopula.mpPopulation;
	mrPopula.mpPopulation = mpNextGen;
	mpNextGen = tmp;

	// Remove the references to the elites from the original
	// population to make the new population their only owner
	for (int i=0; i<mrPopula.mElites; i++)
		// Since we can cut() only with index, not pointer...
		for (int j=0; j<mrPopula.size(); j++)
			// Compare pointers
			if (mpNextGen->getp(j) == &situation.getOrdered(i)) {
				mpNextGen->cut (j);
				// Move the elite holes to the beginning of the population
				if (mpNextGen->getp (i)) {
					mpNextGen->put (mpNextGen->getp (i), j);
					mpNextGen->cut (i);
				}
			}

	// Remove the old population
	// delete mrPopula.mpPopulation;
}

void EAStrategy::addFeaturesTo (Genome& genome) const {
}

void EAStrategy::print (TextOStream& out) {
	out.printf ("Evolving with strategy (e/u[+,]l) = (%d/%d+%d)\n",
				mrPopula.mElites,
				mrPopula.selParams().muFor(mrPopula.size()),
				mrPopula.size());

	out.printf ("Mutation coefficient=%f (binary), %f (int), %f (double rate), "
				"%f (double variance)\n\n",
				mrPopula.mutRate().binaryRate(),
				mrPopula.mutRate().intRate(),
				mrPopula.mutRate().doubleRate(),
				mrPopula.mutRate().doubleVariance());
}

template<class TYPE>
void checkArray (const Array<TYPE>& arr) {
	for (int i=0; i<arr.size(); i++)
		if (arr.getp(i))
			arr.getp(i)->check ();
}

void EAStrategy::check () const {
	if (mpNextGen) {
		mpNextGen->check ();
		checkArray<Individual> (*mpNextGen);
	}
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//       ----- o                              ----                           //
//       |        |    _    ___   ____  ____ (      |   ___   |   ____       //
//       |---  | -+- |/ \  /   ) (     (      ---  -+-  ___| -+- (           //
//       |     |  |  |   | |---   \__   \__      )  |  (   |  |   \__        //
//       |     |   \ |   |  \__  ____) ____) ___/    \  \__|   \ ____)       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void FitnessStats::reset () {
	mMinFitness = 1E30;
	mMaxFitness = 0.0;
	mSumFitness = 0.0;
	mAvgOver = 0;
}

void FitnessStats::add (double fitness) {
	if (fitness < mMinFitness)
		mMinFitness = fitness;

	if (fitness > mMaxFitness)
		mMaxFitness = fitness;

	mSumFitness += fitness;

	mAvgOver++;
}

void FitnessStats::print (TextOStream& out) const {
	out.printf ("Fitness min/avg/max = %f / %f / %f\n",
				mMinFitness, mSumFitness/mAvgOver, mMaxFitness);
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               ----                  |           o                        //
//               |   )       --        |  ___   |           _               //
//               |---   __  |  ) |   | |  ___| -+- |  __  |/ \              //
//               |     /  \ |--  |   | | (   |  |  | /  \ |   |             //
//               |     \__/ |     \__! |  \__|   \ | \__/ |   |             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
* Standard constructor.
*
* @param params["autoAdapt"] Should self-adaptation of mutation rates be used? [Default:0]
* @param params["boolRate"] The mutation rate for @ref BinaryGene genes.
* @param params["intRate"] The mutation rate for @ref IntGene genes.
* @param params["floatRate"] The mutation rate for @ref FloatGene genes.
* @param params["floatVariance"] The mutation variance for @ref FloatGene genes.
*******************************************************************************/
Population::Population (EAEnvironment&   envir,   /**< The environment in which the population evolves in. */
						const StringMap& params)  /**< Additional dynamic parameters as a @ref String @ref Map. */
		: rpEnvironment (&envir)
{
	//
	// Set mutation rates
	//
	double brate = getOrDefault(params,"Population.boolRate", String(0.01)).toDouble ();
	double irate = getOrDefault(params,"Population.intRate", String(0.01)).toDouble ();
	double frate = getOrDefault(params,"Population.floatRate", String(0.1)).toDouble ();
	double fvar =  getOrDefault(params,"Population.floatVariance", String(0.1)).toDouble ();

	mGlobalMutationRate.binaryRate (brate);
	mGlobalMutationRate.intRate (irate);
	mGlobalMutationRate.doubleRate (frate);
	mGlobalMutationRate.doubleVariance (fvar);
	mGlobalMutationRate.autoAdaptation (getOrDefault(params,"Population.autoAdapt", String(0)).toInt ());
	mAutoadjustGMR = false;
}

void Population::check () const
{
}
