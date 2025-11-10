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

#include <magic/mpararr.h>
#include <magic/mdatastream.h>
#include "nhp/simplepopula.h"
#include "nhp/gaenvrnmt.h"
#include "nhp/mutrecord.h"


SimplePopulation::SimplePopulation (EAEnvironment& envir, const StringMap& params)
		: Population (envir, params)
{
	// Set selection parameters
	int popSize = getOrDefault (params, "SimplePopulation.size", String(20)).toInt ();
	mElites     = getOrDefault (params, "EAStrategy.elites", String(int (popSize*0.1))).toInt ();
	mUseGlobalElites = true;

	int u = getOrDefault (params, "Selection.mu", String(int(popSize*0.2))).toInt ();
	if (u>0)
		mSelectionParams.setMu (u);
	else if (mElites>=0)
		mSelectionParams.setMu (u);
	else
		mSelectionParams.setMuPart (0.2);

	// Optional parameter
	mSelectionParams.setMuPart (getOrDefault (params, "Selection.muPart", String(0.2)).toDouble ());

	int    q  = getOrDefault(params,"Selection.q", String(3)).toInt ();
	double ep = getOrDefault(params,"Selection.eta+", String(1.2)).toDouble ();	

	mSelectionParams.setQ (q);
	mSelectionParams.setEtaPlus (ep);

	// Check which strategy parameters to self-adapt
	bool am  = getOrDefault(params,"Selection.adaptMu", String(0)).toInt ();
	bool aq  = getOrDefault(params,"Selection.adaptQ", String(0)).toInt ();
	bool aep = getOrDefault(params,"Selection.adaptEta+", String(0)).toInt ();
	mSelectionParams.adaptParams (am, aq, aep);
	mUseGlobalMu = !am;
	mUseGlobalQ = !aq;
	mUseGlobalEtaPlus = !aep;

	/**************************************************************************/
	// Create a template for an Individual

	Genome templ;
	rpEnvironment->addFeaturesTo (templ);     // Add environment genes
	addFeaturesTo (templ);                    // Add population  genes
	templ.addPrivateGenes (templ, params);    // Add genome      genes
	Individual::addGenesTo (templ, params);   // Add individual  genes
	
	// Set individual-based autoadaptive mutation
	templ.selfadjust (mGlobalMutationRate.autoAdaptation());

	/**************************************************************************/
	// Create the population from the template individual

	// Create population
	mpPopulation = new Array<Individual> ();
	mpPopulation->make (popSize); // 20021125: This was popSize-1 for unknown reason

	// Create individuals
	for (int i=0; i<mpPopulation->size(); i++) {
		mpPopulation->put (new Individual (templ), i);
		(*mpPopulation) [i].setSelector(mSelectionParams);
		(*mpPopulation) [i].init ();
		(*mpPopulation) [i].incarnate (true);
	}

	// Set evolution strategy
	mpStrategy = new EAStrategy (*this);

	failtrace_begin;
	params.failByThrowOnce ();
	if (!params.getp("EAStrategy.silent") || params["EAStrategy.silent"] == "0") {
		mpStrategy->print (mOuts);
	}
	failtrace_end;
	
	// Initialize other parameters;
	minsimilarity = getOrDefault (params, "EAStrategy.minSimilarity", String(0.1)).toDouble ();
	mAge = 0;

	// Set logging
	mEvolog.autoFlush ();
	mOuts.autoFlush ();
}

SimplePopulation::~SimplePopulation () {
	delete mpStrategy;
	delete mpPopulation;
}

void SimplePopulation::addFeaturesTo (Genome& genome) const {
}

double SimplePopulation::evolve (int gens, const char* logfile, double target_fitn)
{
	FUNCTION_BEGIN;
	
	// Open logfile
	FILE* save = NULL;
	if (logfile) {
		save = fopen (logfile, "w");
		ASSERTWITH (save, format ("Log file '%s' couldn't be opened", logfile));

		fprintf (save, "Generation, min_fitness, avg_fitness, max_fitness\n");
		fflush (save);
	}
	if (save)
		mEvolog.setDevice (new File (save));

	// mOuts.setFlag (EAStrategy::TRACE_RECOMBINATION);
	// mOuts.setFlag (EAStrategy::TRACE_MUTATION);
	
	// For a number of generations
	for (int g=0; g<gens; g++) {

		if (MutabilityRecord::record)
			MutabilityRecord::reset ();

		// Evolve for one generation
		failtrace (mpStrategy->evolve (*rpEnvironment, mOuts, mEvolog));

		// Check the termination criteria
		if (target_fitn != -1 && rpEnvironment->bestfitn < target_fitn)
			break;
		
		mEvolog << "\n";
	}
	
	if (save)
		fclose (save);

	FUNCTION_END;
	return mFitnessStats.minFitness();
}

/*******************************************************************************
 * Worker thread to evaluate individuals in a population.
 ******************************************************************************/
class EvaluationWorker : public Thread {
	SimplePopulation*	rpPopula;
	int					mIndividual;
	EAEnvironment*		rpEnvironment;
	TextOStream*		rpOut;
	
  public:
	EvaluationWorker();
	EvaluationWorker(SimplePopulation& popula, int individual, EAEnvironment& environment, TextOStream& out);
	virtual void*		execute		();
};


EvaluationWorker::EvaluationWorker()
		: rpPopula(NULL), mIndividual(0), rpEnvironment(NULL), rpOut(NULL)
{
	fprintf (stderr, "EvaluationWorker::EvaluationWorker() called.\n");
	MUST_OVERLOAD;
}

/*******************************************************************************
 * Initialize evaluation worker by storing callback data.
 ******************************************************************************/
EvaluationWorker::EvaluationWorker(
	SimplePopulation& popula,
	int               individual,
	EAEnvironment&    environment,
	TextOStream&      out)
		: rpPopula(&popula), mIndividual(individual), rpEnvironment(&environment),
		  rpOut(&out)
{
	// out.printf("Created worker %02d.\n", individual);
}

/*******************************************************************************
 * Evaluate the individual.
 ******************************************************************************/
void* EvaluationWorker::execute ()
{
	// fprintf (stderr, "EvaluationWorker::execute () called.\n");
	// rpOut->printf("Evaluating %02d...\n", mIndividual);
	// rpOut->flush();
	rpPopula->evaluate(mIndividual, *rpEnvironment, *rpOut);
	return NULL;
}

void SimplePopulation::evaluate (int i, EAEnvironment& environment, TextOStream& out)
{
	FUNCTION_BEGIN;
		
	double fitness =  (*this) [i].evaluate (environment);

	mThreadLock.lock();
	mFitnessStats.add (fitness);

#if 0
	out.printf ("Indv%3d: ", i);
	failtrace ((*this) [i].print (out));
	out.printf ("", fitness);
	if ((*this)[i].averaged_over ()>1)
		out.printf (" (%d evls)", (*this)[i].averaged_over ());
	if ((*this)[i].getkings())
		out.printf (", kings=%3d\n", (*this)[i].getkings());
	out << "\n";
#endif
	
	mThreadLock.unlock();
	FUNCTION_END;
}

/*******************************************************************************
 *
 ******************************************************************************/
void SimplePopulation::evaluate (EAEnvironment& environment, TextOStream& out)
{
	FUNCTION_BEGIN;

	mFitnessStats.reset ();

	// Evaluate individuals
	// @TODO Evaluate individuals in worker threads
	Array<EvaluationWorker> workers(size());
	for (int i=0; i<size(); i++) {
		EvaluationWorker* worker = new EvaluationWorker(*this, i, environment, out);
		workers.put(worker, i);
		// out.printf("Starting thread %02d\n", i);
		workers[i].start();
	}

	// Wait for threads to end
	for (int i=0; i<size(); i++) {
		// out.printf("Waiting for thread %02d...\n", i);
		workers[i].join();
	}

	out.printf ("SimplePopulation report gen %d: ", mAge);
	mFitnessStats.print (out);

	if (mAutoadjustGMR) {
		// Something here
	}
	
	mAge++;

	FUNCTION_END;
}

void SimplePopulation::resetFitnesses ()
{
	for (int i=0; i<size(); i++) {
	}
}

void SimplePopulation::print (TextOStream& out) const {
	FUNCTION_BEGIN;
	
	out << "SimplePopulation {\n";
	out.printf ("size=%d,\n", size());
	
	for (int i=0; i<size(); i++) {
		(*this)[i].print (out);
		out << "\n";
	}

	out << "}\n";

	FUNCTION_END;
}

void SimplePopulation::report (TextOStream& log) const
{
	log.printf ("%d %.30f %.30f %.30f ", mAge,
				mFitnessStats.minFitness(), mFitnessStats.avgFitness(),
				mFitnessStats.maxFitness());
	if (MutabilityRecord::record)
		log.printf ("%f %f %f %f %f %f %f %.30f %f",
					MutabilityRecord::boolMin(),
					MutabilityRecord::boolAvg(),
					MutabilityRecord::boolMax(),
					MutabilityRecord::floatMin(),
					MutabilityRecord::floatAvg(),
					MutabilityRecord::floatMax(),
					MutabilityRecord::floatVarMin(),
					MutabilityRecord::floatVarAvg(),
					MutabilityRecord::floatVarMax());
	log.flush ();
}

void SimplePopulation::check () const {
	ASSERT (mpPopulation);
	for (int i=0; i<mpPopulation->size(); i++)
		(*mpPopulation)[i].check ();
	rpEnvironment->check ();
	mpStrategy->check ();
}

