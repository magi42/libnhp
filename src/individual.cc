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

#include <magic/mstream.h>
#include <magic/mclass.h>

#include "nhp/gaenvrnmt.h"
#include "nhp/individual.h"
#include "nhp/genes.h"
#include "nhp/population.h"
#include "nhp/selection.h"

impl_dynamic (Individual, {Comparable});

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//               ---           | o       o     |             |               //
//                |    _       |               |        ___  |               //
//                |  |/ \   ---| | |   | |  ---| |   |  ___| |               //
//                |  |   | (   | |  \ /  | (   | |   | (   | |               //
//               _|_ |   |  ---| |   V   |  ---|  \__!  \__| |               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

Individual::Individual () {
	mpSelector = new Selector ();
	incarnate (false);	
}

Individual::Individual (const Individual& prototype) : genome (prototype.genome) {
	mpSelector = new Selector (*prototype.mpSelector);
	incarnate (false);
}

Individual::Individual (const Genome& genotype) : genome (genotype) {
	mpSelector = new Selector ();
	incarnate (false);
}

Individual::~Individual () {
	delete mpSelector;
}

void Individual::addGenesTo (Genome& g, const StringMap& params) {
	Selector::addGenesTo (g, params);
}

void Individual::incarnate (bool doinit) {
	fitness = 0.0;
	avg_over = 0;
	age = 0;

	mFeatures.empty ();
	// mFeatures.set ("age", new Int (0));
	
	if (doinit) {
		// Read selection method parameters from genome
		mpSelector->read (genome);
		
		// Then, launch the ontogenesis of the individual's corpus (it
		// actualizes only if there is the "init" gene in the genome)
		genome.execute (GeneticMsg ("init", *this));
	}
}

void Individual::grow_older () {
	//((long&)static_cast <Int&> (mFeatures ["age"]))++;
	age++;
}


/*******************************************************************************
* Evaluates the fitness of the phenotype in the given environment.
*
* The evaluation is done as many times as defined for the environment. The
* resulting fitness of the evaluations is averaged. The multiple evaluation
* is useful if the environment or the decoding from genotype to phenotype
* has randomness.
*
* The fitness is also cached for later comparisons, if the genotype does
* not change.
*******************************************************************************/
double Individual::evaluate (EAEnvironment& envr, /**< Environment to evaluate the fitness in. */
							 bool           force /**< Force re-evaluation of the fitness.     */) 
{
	// Evaluate as many times as required for averaging.
	while (avg_over < envr.evals() || force) {
		double measured_fitness = envr.evaluate (*this);

		// Take an average of old measurements
 		fitness = (fitness*avg_over + measured_fitness) / (++avg_over);

		if (avg_over >= envr.evals() && force)
			break;
	}

	// Evaluating is tiring
	grow_older ();

	return fitness;
}

void Individual::joinfitness (const Individual& other) {
	if (avg_over>0 || other.avg_over>0)
		fitness = (fitness*avg_over + other.fitness*other.avg_over) / (avg_over+other.avg_over);
	
	avg_over += other.avg_over;

	// Join the number of kinghoods
	for (int k=0; k<other.genome.getkings(); k++)
		genome.addking ();
}


/*
double Individual::weighSelection (const SelectionMatrix& selmat, int i, int j) const {
	return mpSelector->select (selmat, i, j);
}
*/

int Individual::compare (const Comparable& other) const {
	const Individual* o = static_cast <const Individual*> (&other);
	if (fitness < o->fitness)
		return -1;
	else
		if (fitness == o->fitness)
			return 0;
		else
			return 1;
}

void Individual::recombine (const Individual& a, const Individual& b) {
	genome.recombine (a.genome, b.genome);
	incarnate (false);
}

bool Individual::pointMutate (const MutationRate& k) {
	// Forward the request
	bool mut = genome.pointMutate (k);

	// If a mutation has actualized, we are considered a new individual
	if (mut)
		incarnate (false);
	
	return mut;
}

/*
Genstruct* Individual::replicate () const {
	Individual* result = new Individual ();
}
*/

void Individual::print (TextOStream& out) const
{
	FUNCTION_BEGIN;

	if (out[PRINT_CLASSNAMES])
		out << "Indv";
	out << '{';

	if (out[PRINT_NAMES])
		out << "ftn=";
	
	if (fitness>0.0)
		out.printf ("%01.3f, ", fitness);
	else
		out << "[UNKN], ";

	if (out[PRINT_NAMES])
		out << "age=";

	out.printf ("%03d, ", getage ());

	// Print the genome recursively
	failtrace (genome.print (out));

	out << '}';

	FUNCTION_END;
}

DataOStream& Individual::operator>> (DataOStream& out) const {
	out.name ("fitness") << fitness;
	out.name ("avg_over") << avg_over;
	out.name ("age") << age;
	out.name ("genome") << genome;

	return out;
}

void Individual::check () const {
	forMap(String,Object,mFeatures,i) {
		ASSERT (!isempty(i.key()));
		i.value().check ();
	}

	ASSERT (avg_over>=0);
	ASSERT (age>=0);
	genome.check ();
}

void Individual::setSelector (const SelectionPrms& templ) {
	delete mpSelector;
	mpSelector = new Selector (templ);
}
