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

#include "nhp/individual.h"
#include "nhp/population.h"
#include "nhp/simplepopula.h"
#include "nhp/selection.h"
#include "nhp/genes.h"
#include "nhp/mutator.h"
#include <magic/mmath.h>

SelectionMatrix::SelectionMatrix (const SelectionSituation& situation) {
	calculateMatrix (situation);
}

void SelectionMatrix::calculateMatrix (const SelectionSituation& situation) {

	//
	// Compute selection matrices for each method
	//

	// Init the final matrix
	int popSize = situation.population().size();
	mSelection.make (popSize, popSize);
	mSelection = 0.0;
	
	double rowsum;
	for (int i=0; i<popSize; i++) {
		rowsum=0.0;
		
		// Let the individual tell it's selection willingness towards
		// every individual, including itself
		for (int j=0; j<popSize; j++)
			rowsum += (mSelection.get(i,j) = situation.getOrdered(i).selector().select(situation, i, j));
		
		// Equalize row sum to 1.0
		// rowsum = 1/(rowsum*mOPop.size);
		for (int j=0; j<popSize; j++)
			mSelection.get(i,j) /= rowsum;
	}
	// sout << mSelection;
	
	// Sum the matrices
	mSelection.multiplyToSum (1.0);
	mSelection *= transpose (mSelection);
	mSelection.multiplyToSum (1.0);
}

void SelectionMatrix::selectRandomPair (int& a, int& b) const {
	double rn=frnd();
	double pos=0.0;
	for (int i=0; i<mSelection.rows; i++) {
		for (int j=0; j<mSelection.cols; j++) {
			if ((pos+=mSelection.get(i,j)) >= rn) {
				// TRACE2 ("Selected %d+%d", i,j);
				a = i;
				b = j;
				return;
			}
		}
	}
	FORBIDDEN;
}



/////////////////////////////////////////////////////////////////////////////////////////
//  ----       |                 o             ---- o                     o            //
// (      ___  |  ___   ___   |           _   (        |         ___   |           _   //
//  ---  /   ) | /   ) |   \ -+- |  __  |/ \   ---  | -+- |   |  ___| -+- |  __  |/ \  //
//     ) |---  | |---  |      |  | /  \ |   |     ) |  |  |   | (   |  |  | /  \ |   | //
// ___/   \__  |  \__   \__/   \ | \__/ |   | ___/  |   \  \__!  \__|   \ | \__/ |   | //
/////////////////////////////////////////////////////////////////////////////////////////

SelectionSituation::SelectionSituation (const SimplePopulation& pop)
		: mrPop (pop), mOrdPop (pop.getPopArray()) {
	mOrdPop.quicksort ();
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ----       |                 o            ----                      //
//     (      ___  |  ___   ___   |           _   |   )            ____     //
//      ---  /   ) | /   ) |   \ -+- |  __  |/ \  |---  |/\ |/|/| (         //
//         ) |---  | |---  |      |  | /  \ |   | |     |   | | |  \__      //
//     ___/   \__  |  \__   \__/   \ | \__/ |   | |     |   | | | ____)     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

SelectionPrms::SelectionPrms () {
	mMu = -1;
	mMuPart = 0.2;
	mEtaPlus = 1.0;
	mQ = 3;
	mAdaptiveMu = false;
	mAdaptiveEtaPlus = false;
	mAdaptiveQ = false;

	mWeightedSelection = false;

	mSelMethodW.make (Selector::number_of_methods);
	useOnlyMethod (Selector::MULAMBDASELECTION);
}

void SelectionPrms::adaptParams (bool adapt_mu, bool adapt_etaplus, bool adapt_q) {
	ASSERTWITH (!adapt_mu || mMu<0,
				"Must use partial mu with self-adaptation");
	
	mAdaptiveMu = adapt_mu;
	mAdaptiveEtaPlus = adapt_etaplus;
	mAdaptiveQ = adapt_q;
}

void SelectionPrms::useOnlyMethod (int m) {
	ASSERT (m>=0 && m<Selector::number_of_methods);

	for (int i=0; i<Selector::number_of_methods; i++)
		mSelMethodW[i] = 0.0;

	mSelMethodW[m] = 1.0;

	mAdaptiveWeights = false;
}

void SelectionPrms::setEtaPlus (double etaplus) {
	ASSERT (etaplus>=1 && etaplus<=2);
	mEtaPlus = etaplus;
}

void SelectionPrms::setQ (int q) {
	ASSERT (q>1);
	mQ = q;
}

void SelectionPrms::setMu (int m) {
	ASSERT (m>0);
	mMu = m;
	mMuPart = -1;
	mAdaptiveMu = false;
}

void SelectionPrms::setMuPart (double muPart) {
	ASSERT (muPart>=0 && muPart<1);
	mMu = -1;
	mMuPart = muPart;
}

int SelectionPrms::muFor (int populsize) const {
	if (mMu>=0)
		return mMu;
	else {
		ASSERT (mMuPart>=0);
		return int (mMuPart*populsize+0.5);
	}
}

void SelectionPrms::copy (const SelectionPrms& o) {
	mMu                = o.mMu;
	mMuPart            = o.mMuPart;
	mEtaPlus           = o.mEtaPlus;
	mQ                 = o.mQ;
	mAdaptiveMu        = o.mAdaptiveMu;
	mAdaptiveEtaPlus   = o.mAdaptiveEtaPlus;
	mAdaptiveQ         = o.mAdaptiveQ;

	mWeightedSelection = o.mWeightedSelection;
	mSelMethodW        = o.mSelMethodW;
	mAdaptiveWeights   = o.mAdaptiveWeights;
}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                   ----       |                                           //
//                  (      ___  |  ___   ___   |                            //
//                   ---  /   ) | /   ) |   \ -+-  __  |/\                  //
//                      ) |---  | |---  |      |  /  \ |                    //
//                  ___/   \__  |  \__   \__/   \ \__/ |                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

Selector::Selector (const SelectionPrms& orig) : SelectionPrms (orig), mScore (0) {
}

// We must use our own mutator that is not affected by the possible
// self-adaptation of the mutation rates
ConstantMutator selectorMutator (0.02);

void Selector::addGenesTo (Gentainer& g, const StringMap& params) {

	// Parameters of different selection methods

	// mMuPart (for u% selection)
	//g.add (&(new FloatGene ("u%", 0.01, 0.99, 0.2))->setMutator(&selectorMutator));
	g.add (new BitFloatGene ("u%", 0.01, 0.99, 16, params));
	g.add (new IntGene ("q", 2, 10, 1.0));		// mQ (for tournament selection)
	g.add (new FloatGene ("e+", 1, 2, 0.2));	// mEtaPlus (for rank-based selection)

	g["u%"].hide ();
	g["q"].hide ();
	g["e+"].hide ();

	// Weights/propabilities of the different selection methods
	for (int i=0; i<4; i++) {
		g.add (new FloatGene (format ("s%d",i), 0, 1, 0.2));
		g[(CONSTR)format ("s%d",i)].hide ();
	}
}

void Selector::read (const Genome& g) {
	try {
		// Selection method parameters
		if (mAdaptiveMu)
			mMuPart = static_cast <const FloatGene*> (g.getGene ("u%"))->getvalue();
		if (mAdaptiveEtaPlus)
			mEtaPlus = static_cast <const FloatGene*> (g.getGene ("e+"))->getvalue();
		if (mAdaptiveQ)
			mQ = static_cast <const IntGene*> (g.getGene ("q"))->getvalue();
		
		// Selection method weights/propabilities
		if (mAdaptiveWeights) {
			mSelMethodW.make (Selector::number_of_methods);
			for (int i=0; i<Selector::number_of_methods; i++) {
				const Genstruct& gene = *g.getGene (format ("s%d", i));
				mSelMethodW[i] = dynamic_cast <const FloatGene&> (gene).getvalue();
			}
			// Ensure that selection method weights sum to 1.0
			multiplyToUnity (mSelMethodW);
		}
	} catch (...) {
		ASSERT (false);
	}
}

double Selector::select (const SelectionSituation& situation, int i, int j) const {
	if (mWeightedSelection) {
		// Add the results of different selection functions together,
		// weighing them nicely
		double love=0.0;
		for (int m=0; m<Selector::number_of_methods; m++)
			love += mSelMethodW[m] * selectWithMethod (situation, m, i, j);
		return love;
	} else {
		//
		// Choose one selection method using the propabilities for
		// different methods
		//
		double p=frnd ();
		for (int m=0; m<Selector::number_of_methods; p-=mSelMethodW[m++])
			if (p <= mSelMethodW[m])
				return selectWithMethod (situation, m, i, j);
		FORBIDDEN; // Probibility selection is not supported currently
	}
}

double Selector::selectWithMethod (const SelectionSituation& situation, int sm, int i, int j) const {
	ASSERTWITH (sm>=0 && sm<=3, "Selection method index out of range");

	switch (sm) {
	  case 0: return muLambdaSelection (situation, i, j);
	  case 1: return linearRanking (situation, i, j);
	  case 2: return proportionalSelection (situation, i, j);
	  case 3: return tournamentSelection (situation, i, j);
	};

	/* TODO: Handle error situation. */
	return 0.0;
}

double Selector::muLambdaSelection (const SelectionSituation& situation, int i, int j) const {
	const SimplePopulation& pop = situation.population();
	int mu;
	if (pop.mUseGlobalMu)
		mu = pop.selParams().muFor (pop.size());
	else // Use locally autoadapted mu
		mu = muFor (pop.size());

	return (j>mu)? 0.0 : 1.0;
}

double Selector::linearRanking (const SelectionSituation& situation, int i, int j) const {
	double ep = situation.population().mUseGlobalEtaPlus?
		situation.population().selParams().mEtaPlus
		: mEtaPlus;

	return (ep-(2*ep-2)*double(j-1)/double(situation.population().size()-1));
}

double Selector::proportionalSelection (const SelectionSituation& situation, int i, int j) const {
	return situation.getOrdered(j).getfitness() / situation.population().mFitnessStats.avgFitness();
}

double Selector::tournamentSelection (const SelectionSituation& situation, int i, int j) const {
	double q = situation.population().mUseGlobalQ?
		int(situation.population().selParams().mQ)
		: mQ;
	int popSize = situation.population().size();
	return pow (popSize, -q) * (pow(popSize-j+1, q)-pow(popSize-j, q)); 
}


/*
void Selector::operator= (const Selector& o) {
	SelectionPrms::operator= (o);
	mScore = o.mScore;
}

void Selector::operator= (const SelectionPrms& o) {
	SelectionPrms::operator= (o);
	mScore = 0;
}
*/
