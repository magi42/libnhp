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
#include <magic/mclass.h>
#include <magic/mtextstream.h>

#include "nhp/genetics.h"
#include "nhp/genes.h"
#include "nhp/mutrecord.h"
#include "nhp/mutator.h"

impl_dynamic (Genstruct, {Object});
impl_dynamic (Gentainer, {Genstruct});
impl_dynamic (Genome, {Gentainer});



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//       |   |                     o            ----                        //
//       |\ /|        |   ___   |           _   |   )  ___   |   ___        //
//       | V | |   | -+-  ___| -+- |  __  |/ \  |---   ___| -+- /   )       //
//       | | | |   |  |  (   |  |  | /  \ |   | | \   (   |  |  |---        //
//       |   |  \__!   \  \__|   \ | \__/ |   | |  \   \__|   \  \__        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

MutationRate::MutationRate (const Gentainer& g) {
	mBinaryRate	= static_cast<const FloatGene&> (*g.getGene ("Rb")).getvalue();
	mFloatRate		= 1;//static_cast<FloatGene&> (*g.getGene ("Rf")).getvalue();
	mFloatVariance	= static_cast<const FloatGene&> (*g.getGene ("Vf")).getvalue();
	mIntRate		= static_cast<const FloatGene&> (*g.getGene ("Ri")).getvalue();
	mOneBitMutation = false;
}

MutationSimpleMutator mutmut (30);

void MutationRate::addGenesTo (Gentainer& g, const StringMap& params)
{
	double low = isnull(params)? 0.01 : getOrDefault (params, "MutationRate.lowBound", String(0.01)).toDouble ();
	g.add (&(new FloatGene ("Rb", low, 1, 10.0))->setMutator(&mutmut));
	g.add (&(new FloatGene ("Rf", low, 1, 10.0))->setMutator(&mutmut));
	//g.add (&(new FloatGene ("Vf", 0, 1, 10.0))->setMutator(&mutmut));
	g.add (&(new FloatGene ("Vf", 0, 1, 10.0))->setMutator(&mutmut));
	//g.add (new FloatGene ("Vf", 0, 1, 10.0));
	g.add (&(new FloatGene ("Ri", low, 1, 10.0))->setMutator(&mutmut));

	// Set them hidden by default
	g["Rb"].hide ();
	g["Rf"].hide ();
	g["Vf"].hide ();
	g["Ri"].hide ();
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//               ----                                                        //
//              |      ___    _    ____  |             ___   |               //
//              | --- /   ) |/ \  (     -+- |/\ |   | |   \ -+-              //
//              |   \ |---  |   |  \__   |  |   |   | |      |               //
//              |___/  \__  |   | ____)   \ |    \__!  \__/   \              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

Genstruct::Genstruct (const GeneticID& nam) {
	size = 0;
	id = nam;
	mHidden = false;
}

Genstruct::Genstruct (const Genstruct& orig) {
	size = orig.size;
	id = orig.id;
	mHidden = orig.mHidden;
}

int Genstruct::length () const {
	if (!size)
		size = calc_len ();
	
	return size;
}

void Genstruct::print (TextOStream& out) const {
	if (!mHidden)
		out.printf ("%s", (CONSTR) id);
}

DataOStream& Genstruct::operator>> (DataOStream& out) const {
	out.name ("size") << size;
	out.name ("id") << id;
	out.name ("mHidden") << int(mHidden);
	return out;
}

void Genstruct::check () const {
	ASSERT (size>=0);
	ASSERTWITH (size<1000000, "A sensible upper limit for error checking");
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                ----                       o                               //
//               |      ___    _    |   ___      _    ___                    //
//               | --- /   ) |/ \  -+-  ___| | |/ \  /   ) |/\               //
//               |   \ |---  |   |  |  (   | | |   | |---  |                 //
//               |___/  \__  |   |   \  \__| | |   |  \__  |                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

Gentainer::Gentainer (const GeneticID& iid) : Genstruct (iid) {
	self_adjust = false;
	mRecombRate = 1.0;
}

Gentainer::Gentainer (const Gentainer& orig) : Genstruct (orig) {
	// Replicate all the other's substructures
	for (int i=0; i<orig.substructs.size(); i++)
		substructs.add (orig.substructs[i].replicate ());
	
	self_adjust = orig.self_adjust;
	mRecombRate = orig.mRecombRate;
}

void Gentainer::add (Genstruct* genestr) {
	ASSERTWITH (genestr, "Genstruct to be added to Gentainer must not be null pointer.");
	substructs.add (genestr);
}

void Gentainer::init () {
	// Spread the initialization message to all substructures
	for (int i=0; i<substructs.size(); i++)
		substructs[i].init ();
}

void Gentainer::addPrivateGenes (Gentainer& parent, const StringMap& params) {
	// Propability of crossovers
	if (isnull(params))
		add (new FloatGene ("Px", 0, 1, 0.1)); // Same as below!
	else  {
		params.failByNullOnce ();
		if (isnull (params["Gentainer.recombFreq"]))
			add (new FloatGene ("Px", 0, 1, 0.1)); // Same as above!
		else {
			params.failByNullOnce ();
			double value = params["Gentainer.recombFreq"].toDouble ();
			add (new FloatGene ("Px", value, value, 0.0));
		}
	}
	
	add (new IntGene ("RM", 1, 2, 1.0));	// Recombination method (1=crossover, 2=uniform)
	add (new IntGene ("Nx", 1, 4, 1.0));	// Number of crossovers

	// Hide 'em
	//(*this)["RM"].hide ();
	//(*this)["Nx"].hide ();
	//(*this)["Px"].hide ();

	MutationRate::addGenesTo (*this, params);

	// Also add private genes for contained genes
	for (int i=0; i<substructs.size(); i++)
		substructs[i].addPrivateGenes (*this, params);

	// (Well, most of them don't propably exist at this time yet)
}

const Genstruct* Gentainer::getGene (const GeneticID& cc_name) const
{
	try {
		// First seek this gentainer width-first
		for (int i=0; i<substructs.size(); i++) {
			if (substructs[i].getID() == cc_name)
				return &substructs[i];
		}

		// The gene might be somewhere below, so now we seek recursively
		for (int i=0; i<substructs.size(); i++)
			if (const Genstruct* gotten = substructs[i].getGene (cc_name))
				return gotten;
		
	} catch (exception e) {
		// Forward any exceptions
		throw invalid_gene_name (format ("Error while seeking gene '%s':\n%s",
										 (CONSTR)cc_name, e.what()));
	}

	// Gene really not found
	return NULL;
}

bool Gentainer::pointMutate (const MutationRate& k) {
	bool mutated = false;

	if (self_adjust) {
		MutationRate mult (*this);
		MutationRate combined (k, mult);

		MutabilityRecord::addBoolMutability (combined.binaryRate());
		MutabilityRecord::addFloatMutability (combined.doubleRate());
		MutabilityRecord::addFloatVariance (combined.doubleVariance());
		
		for (int i=0; i<substructs.size(); i++)
			if (substructs[i].pointMutate (combined))
				mutated = true;
	} else {
		// No self-adjustment
		for (int i=0; i<substructs.size(); i++)
			if (substructs[i].pointMutate (k))
				mutated = true;
	}

	return mutated;
}

void Gentainer::recombine (const Genstruct& as, const Genstruct& bs) {
	const Gentainer& a = static_cast<const Gentainer&> (as);
	const Gentainer& b = static_cast<const Gentainer&> (bs);

	// Currently the genomes have to have equal length
	ASSERTWITH (a.length() == length() && b.length() == length(),
				format ("Parameter error, len=%d, a.len=%d, b.len=%d",
						length(), a.length(), b.length()));
	
	// Make an array of cross-over position marks
	bool crosspos [substructs.size()-1];
	for (int i=0; i<substructs.size()-1; i++)
		crosspos [i] = false;
	
	// Mark the cross-over points
	int n = static_cast<const AnyIntGene&> (*getGene ("Nx")).getvalue();
	double pX = static_cast<const AnyFloatGene&> (*getGene ("Px")).getvalue();
	for (int i=0; i<n; i++)
		if (frnd ()<pX)
			crosspos [rnd (substructs.size())] = true;

	// And copy the parents, switching at the marks...
	int whichpar = 0;
	for (int i=0; i<substructs.size(); i++) {
		ASSERTWITH (substructs[i].getID() == a[i].getID() &&
					substructs[i].getID() == b[i].getID(),
					"Chromosomes must be in equal order to be crossed");

		// If at crossover mark
		if (i>0 && crosspos [i-1]) {
			// Switch order
			whichpar = 1-whichpar;
			
			// recurse the crossover
			substructs[i].recombine (whichpar? a[i]:b[i], (1-whichpar)? a[i]:b[i]);
		} else
			// otherwise just copy as is
			substructs[i].copy (whichpar? a[i] : b[i]);
	}
}

double Gentainer::equality (const Genstruct& o) const {
	const Gentainer& other = static_cast<const Gentainer&> (o);

	ASSERT (length() == other.length());

	// For each gene
	double tot=0.0;
	for (int i=0; i<substructs.size(); i++) {
		ASSERTWITH (substructs[i].getID() == other[i].getID(),
					"Chromosomes must be in equal order to be compared");

		// Let gene compare itself to the other
		tot += substructs[i].equality (other.substructs[i]);
	}

	// Return difference
	return tot;
}

Genstruct* Gentainer::replicate	() const {
	// Create a gentainer of the same class. Using this scheme we
	// don't have to implement a replication operation for all
	// different gentainers. Hmm. This might be a wrong approach.
	Gentainer* cln = static_cast<Gentainer*> (dyncreate (this->getclassname())); 

	// There's definitely something missing right here
	
	// Replicate all substructures
	for (int i=0; i<substructs.size(); i++)
		cln->add (substructs[i].replicate ());

	return cln;
}

void Gentainer::copy (const Genstruct& o) {
	copyGenstr (o);
	const Gentainer& other = static_cast<const Gentainer&> (o);
	
	ASSERTWITH (length()==o.length() || length()==0,
				format ("Cannot copy structure of differing dimension (this=%d, other=%d)"
						, length(), o.length()));
	ASSERTWITH (substructs.size()==other.substructs.size() || substructs.size()==0,
				format ("Cannot copy structure of differing dimension (this=%d, other=%d)"
						, substructs.size(), other.substructs.size()));

	if (substructs.size()>0) {
		// Copy all substructures
		for (int i=0; i<substructs.size(); i++)
			substructs[i].copy (other.substructs[i]);
	} else {
		// Can't copy, so clone
		for (int i=0; i<other.substructs.size(); i++)
			substructs.add (other.substructs[i].replicate());
	}

	self_adjust = other.self_adjust;
}

void Gentainer::print (TextOStream& out) const
{
	if (!isempty(id))
		out.printf ("%s=", (CONSTR) id);
	out.printf ("%s {", (CONSTR) getclassname());

	// For all contained genstructs
	for (int i=0; i<substructs.size(); i++) {
		// Print genstruct, if visible
		if (!substructs[i].isHidden()) {
			substructs[i].print (out);
			
			// Print whitespace, but not at end of list, and not for bits.
			if (i<substructs.size()-1 && !substructs[i].is_a("BinaryGene") && !substructs[i+1].is_a("BinaryGene"))
				out << ' ';
		}
	}
	out.flush ();
	out << '}';
}

int Gentainer::calc_len () const {
	int sum=0;
	for (int i=0; i<substructs.size(); i++)
		sum += substructs[i].length ();
	return sum;
}

/*******************************************************************************
* Executes a genetic message in the gentainer.
*
* Sends the message to every contained genstruct.
*
* Implementation for @ref Genstruct.
*******************************************************************************/
bool Gentainer::execute (const GeneticMsg& msg) const
{
	FUNCTION_BEGIN;
	
	// Send the message to the named receiver
	for (int i=0; i<substructs.size(); i++)
		if (substructs[i].getID() == msg.mrReceiver)
			return substructs[i].execute (msg);

	// KLUDGE
	/*
	if (msg.receiver != "init")
		throw exception (format ("Receiver '%s' not found in gentainer %s '%s'",
								 (CONSTR) msg.receiver, (CONSTR) getclassname(),
								 (CONSTR) id));
								 */
	if (msg.mrReceiver == "init")
		return true;

	return false;
	FUNCTION_END;
}

void Gentainer::check () const {
	Genstruct::check ();
	ASSERT (mRecombRate>=0 && mRecombRate<=1);
	substructs.check();
}

DataOStream& Gentainer::operator>> (DataOStream& out) const {
	Genstruct::operator>> (out);
	out.name ("substructs") << substructs;
	return out;
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      ----                                                //
//                     |      ___    _               ___                    //
//                     | --- /   ) |/ \   __  |/|/| /   )                   //
//                     |   \ |---  |   | /  \ | | | |---                    //
//                     |___/  \__  |   | \__/ | | |  \__                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

Genome::Genome () {
	kings = 0;
}

Genome::Genome (const Genome& other) : Gentainer (other) {
	kings = other.kings;
}

void Genome::init () {
	kings = 0;
	Gentainer::init ();
}

void Genome::addPrivateGenes (Gentainer& parent, const StringMap& pars) {
	Gentainer::addPrivateGenes (*this, pars);
}

void Genome::print (TextOStream& out) const {
	Gentainer::print (out);
}

DataOStream& Genome::operator>> (DataOStream& out) const {
	Gentainer::operator>> (out);
	out.name ("kings") << kings;
	return out;
}

void Genome::check () const {
	Gentainer::check ();
	ASSERT (kings>=0);
	ASSERTWITH (kings<100000, "Reasonable upper limit");
}

/*
void Genome::recombine (const Genstruct& as, const Genstruct& bs) {
	Gentainer::recombine (as, bs);
}
*/


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ___   o      |      o     |  ----                                   //
//      |  \     --  |            | |      ___    _               ___       //
//      |   | | |  ) |  __  |  ---| | --- /   ) |/ \   __  |/|/| /   )      //
//      |   | | |--  | /  \ | (   | |   \ |---  |   | /  \ | | | |---       //
//      |__/  | |    | \__/ |  ---| |___/  \__  |   | \__/ | | |  \__       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

