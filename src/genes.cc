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

#include "nhp/genes.h"
#include "nhp/individual.h"
#include "nhp/mutrecord.h"
#include "nhp/mutator.h"

impl_dynamic (Gene, {Genstruct});
impl_dynamic (BinaryGene, {Gene});
impl_dynamic (AnyFloatGene, {Gene});
impl_dynamic (FloatGene, {AnyFloatGene});
impl_dynamic (BitFloatGene, {AnyFloatGene});
impl_dynamic (AnyIntGene, {Gene});
impl_dynamic (IntGene, {AnyIntGene});
impl_dynamic (BitIntGene, {AnyIntGene});
impl_dynamic (InterGene, {Gene});



double mutateFloat (double x, double rate, int bits)
{
	// Encode value in bits
	PackArray<int> binary (bits);

	// TODO
	return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           ----                                            //
//                          |      ___    _    ___                           //
//                          | --- /   ) |/ \  /   )                          //
//                          |   \ |---  |   | |---                           //
//                          |___/  \__  |   |  \__                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

Gene::Gene (const GeneticID& n, double mut) : Genstruct (n) {
	mutability = mut;
}

void Gene::copy (const Genstruct& o) {
	copyGenstr (o);
	mutability = static_cast<const Gene&>(o).mutability;
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           ----  o                        ----                             //
//           |   )     _    ___            |      ___    _    ___            //
//           |---  | |/ \   ___| |/\ \   | | --- /   ) |/ \  /   )           //
//           |   ) | |   | (   | |    \  | |   \ |---  |   | |---            //
//           |___  | |   |  \__| |     \_/ |___/  \__  |   |  \__            //
//                                    \_/                                    //
///////////////////////////////////////////////////////////////////////////////

BinaryGene::BinaryGene (const GeneticID& id, double mut, double initP) : Gene (id, mut) {
	ASSERT (initP>=0 && initP<=1);

	mInitP = initP;
	init ();
}

void BinaryGene::init () {
	mValue = (frnd()<mInitP)? 1:0;
}

bool BinaryGene::pointMutate (const MutationRate& mut_rate) {
	const double phi = 0.22; // Hmm, should be defined somewhere else... Oh well...
	
	// Mutate mutability
	if (false && mut_rate.autoAdaptation()) {
		double p=mutability*mut_rate.binaryRate();
		p = 1/(1+(1-p)/p*exp(-phi*gaussrnd(1)));
		if (p<0.01)
			p = 0.01;
		mutability = p/mut_rate.binaryRate();
	}

	//if (MutabilityRecord::record)
	//	MutabilityRecord::addBoolMutability (mutability);
	
	// Mutate value
	if (frnd()>mut_rate.binaryRate()*mutability)
		return false;

	mValue = mValue? 0:1;
	return true;
}

void BinaryGene::copy (const Genstruct& o) {
	Gene::copy (o);
	shallowCopy (static_cast<const BinaryGene&>(o));
}

void BinaryGene::print (TextOStream& out) const {
	out.printf ("%c", mValue? '1':'0');
}

DataOStream& BinaryGene::operator>> (DataOStream& out) const {
	out.name("id") << id;
	out.name("mutability") << mutability;
	out.name("value") << int(mValue);
	out.name("initP") << mInitP;

	return out;
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//       _               ----- |                 ----                       //
//      / \    _         |     |       ___   |  |      ___    _    ___      //
//     /   \ |/ \  \   | |---  |  __   ___| -+- | --- /   ) |/ \  /   )     //
//     |---| |   |  \  | |     | /  \ (   |  |  |   \ |---  |   | |---      //
//     |   | |   |   \_/ |     | \__/  \__|   \ |___/  \__  |   |  \__      //
//                  \_/                                                     //
//////////////////////////////////////////////////////////////////////////////

AnyFloatGene::AnyFloatGene (const GeneticID& id, double mi, double ma, double mut) : Gene (id, mut) {
	ASSERTWITH (mi<=ma, format ("min (was %f) value must be <= than max value (was %f)"
							   , mi, ma));
	mMin = mi;
	mMax = ma;
}

void AnyFloatGene::copy (const Genstruct& o) {
	Gene::copy (o);
	shallowCopy (static_cast<const AnyFloatGene&>(o));
}

void AnyFloatGene::shallowCopy (const AnyFloatGene& o) {
	mMin	= o.mMin;
	mMax	= o.mMax;
}

void AnyFloatGene::check () const {
	ASSERT (mMin<mMax);
	Gene::check ();
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               ----- |                 ----                               //
//               |     |       ___   |  |      ___    _    ___              //
//               |---  |  __   ___| -+- | --- /   ) |/ \  /   )             //
//               |     | /  \ (   |  |  |   \ |---  |   | |---              //
//               |     | \__/  \__|   \ |___/  \__  |   |  \__              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

FloatGene::FloatGene (const GeneticID& id, double mi, double ma, double mut, double val) : AnyFloatGene (id, mi, ma, mut) {
	value = val;
	mVariance = 1;
	mMutator = NULL;
	mCircular = false;

	if (val<mMin || val>mMax)
		init ();
}

void FloatGene::init () {
	value = mMin + frnd ()*(mMax-mMin);
}

void FloatGene::copy (const Genstruct& o) {
	AnyFloatGene::copy (o);
	shallowCopy (static_cast<const FloatGene&>(o));
}

void FloatGene::shallowCopy (const FloatGene& o) {
	value		= o.value;
	mVariance	= o.mVariance;
	mMutator	= o.mMutator;
}

bool FloatGene::pointMutate (const MutationRate& mut_rate) {
	// Mutate mutability
	if (false && mut_rate.autoAdaptation() && frnd()<mut_rate.doubleRate()) {
		mutability = fabs (mutability + gaussrnd (mutability*0.5));
	}

	if (mMutator) {
		value = mMutator->mutate (value, mMin, mMax, mVariance*mut_rate.doubleVariance());
	} else
		if (mMax>mMin) { // Can be 0 -> gene is immutable
			// Default mutation
			if (frnd()<mut_rate.doubleRate()) {
				double delta = gaussrnd (mut_rate.doubleVariance());
				if (value+delta<mMin)
					delta = mMin;
				else if (value+delta>mMax)
					delta = mMax;
				
				value += delta;
			}
		}

	//if (MutabilityRecord::record)
	//	MutabilityRecord::addFloatMutability (mutability);
	
	return true;
}

double FloatGene::equality (const Genstruct& o) const {
	const FloatGene& other = static_cast<const FloatGene&>(o);
	return fabs(value-((FloatGene&)other).value)/(mMax-mMin);
}

void FloatGene::print (TextOStream& out) const {
	out.printf ("%s=%0.2f", (CONSTR) id, value);
}

void FloatGene::check () const {
	AnyFloatGene::check ();
	ASSERT (value>=mMin);
	ASSERT (value<=mMax);
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ----  o     ----- |                 ----                          //
//        |   )    |  |     |       ___   |  |      ___    _    ___         //
//        |---  | -+- |---  |  __   ___| -+- | --- /   ) |/ \  /   )        //
//        |   ) |  |  |     | /  \ (   |  |  |   \ |---  |   | |---         //
//        |___  |   \ |     | \__/  \__|   \ |___/  \__  |   |  \__         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

BitFloatGene::BitFloatGene (const GeneticID& id, double mi, double ma, int bits, const StringMap& params, double mut) : AnyFloatGene (id, mi, ma, mut) {
	ASSERT (bits>=1 && bits<=32);

	mBitCount = bits;
	for (int i=0; i<bits; i++)
		mBits.add (new BinaryGene (format("b%d", i)));

	Gentainer* dummy = NULL;
	mBits.addPrivateGenes (*dummy, params);
	mBits.selfadjust (false);

	trywith (mGrayCoded = getOrDefault (params, "BitFloatGene.graycoding", String(1)).toInt(),
			 format ("Gene name='%s'", (CONSTR) id));
}

void BitFloatGene::copy (const Genstruct& o) {
	AnyFloatGene::copy (o);
	shallowCopy (static_cast<const BitFloatGene&>(o));
}

void BitFloatGene::shallowCopy (const BitFloatGene& o) {
	mBitCount	= o.mBitCount;
	mBits.copy (o.mBits);
	mGrayCoded = o.mGrayCoded;
}

double BitFloatGene::getvalue () const
{
	// Read the bits to a binary vector
	PackArray<int> bits (mBitCount);
	for (int i=0; i<mBitCount; i++)
		bits[i] = static_cast<const BinaryGene&> (*mBits.getGene(format("b%d",i))).getvalue();

	if (mGrayCoded) {
		// Convert from Gray to binary
		for (int i=mBitCount-2; i>=0; i--)
			bits[i] = bits[i+1]^bits[i];
	}

	// Convert from binary to double
	unsigned long int sum = 0;
	for (unsigned long int i=0, m=1; i < (uint) mBitCount; i++, m<<=1)
		sum |= m*bits[i];

	return mMin+(mMax-mMin)*double(sum)/double((1<<mBitCount));
}

void BitFloatGene::print (TextOStream& out) const {
	out.printf ("%s=", (CONSTR) id);
	mBits.print (out);
}

void BitFloatGene::recombine (const Genstruct& ar, const Genstruct& br) {
	const BitFloatGene& a = static_cast<const BitFloatGene&> (ar);
	const BitFloatGene& b = static_cast<const BitFloatGene&> (br);
	
	mBits.recombine (a.mBits, b.mBits);
}

void BitFloatGene::check () const {
	AnyFloatGene::check ();
	ASSERT (getvalue()>=mMin);
	ASSERT (getvalue()<=mMax);
	mBits.check ();
	ASSERT (mBitCount>=1);
	ASSERT (mBitCount<1000); // Reasonable maximum
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//            _               ---            ----                            //
//           / \    _          |    _    |  |      ___    _    ___           //
//          /   \ |/ \  \   |  |  |/ \  -+- | --- /   ) |/ \  /   )          //
//          |---| |   |  \  |  |  |   |  |  |   \ |---  |   | |---           //
//          |   | |   |   \_/ _|_ |   |   \ |___/  \__  |   |  \__           //
//                       \_/                                                 //
///////////////////////////////////////////////////////////////////////////////

AnyIntGene::AnyIntGene (const GeneticID& id, int mi, int ma, double mut) : Gene (id, mut) {
	ASSERTWITH (mi<=ma, format ("min (was %f) value must be <= than max value (was %f)"
							   , mi, ma));
	mMin = mi;
	mMax = ma;
}

void AnyIntGene::copy (const Genstruct& o) {
	Gene::copy (o);
	shallowCopy (static_cast<const AnyIntGene&>(o));
}

void AnyIntGene::shallowCopy (const AnyIntGene& o) {
	mMin	= o.mMin;
	mMax	= o.mMax;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                    ---            ----                                    //
//                     |    _    |  |      ___    _    ___                   //
//                     |  |/ \  -+- | --- /   ) |/ \  /   )                  //
//                     |  |   |  |  |   \ |---  |   | |---                   //
//                    _|_ |   |   \ |___/  \__  |   |  \__                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

IntGene::IntGene (const GeneticID& nam, int mi, int ma, double mut, int val) : AnyIntGene (nam, mi, ma, mut) {
	ASSERTWITH (mi<ma, format ("min (was %f) value must be smaller than max value (was %f)"
							   , mi, ma));
	mValue = val;
	if (val<mMin || val>mMax)
		init ();
}

int IntGene::getvalue () const {
	return mValue;
}

void IntGene::init () {
	mValue = mMin + rnd (mMax-mMin);
}

void IntGene::copy (const Genstruct& o) {
	AnyIntGene::copy (o);
	//copyGenstr (o);
	shallowCopy (static_cast<const IntGene&>(o));
}

void IntGene::shallowCopy (const IntGene& o) {
	mValue		= o.mValue;
	mMin		= o.mMin;
	mMax		= o.mMax;
}

bool IntGene::pointMutate (const MutationRate& mut_rate) {
	switch (0) {
	  case 0: {
		  if (frnd()<=mut_rate.intRate()*mutability)
			  init ();
	  } break;
	  case 1: {
		  int delta;
		  do {
			  delta = int (gaussrnd (mut_rate.doubleVariance()*mutability));
		  } while (mValue+delta<mMin || mValue+delta>mMax);
		  mValue += delta;
	  }
	};
	return true;
}

double IntGene::equality (const Genstruct& o) const {
	const IntGene& other = static_cast<const IntGene&>(o);
	return double(abs(mValue-((IntGene&)other).mValue))/double(mMax-mMin);
}

void IntGene::print (TextOStream& out) const {
	out.printf ("%s=%d", (CONSTR) id, mValue);
}

void IntGene::check () const {
	AnyIntGene::check ();
	ASSERT (mValue>=mMin);
	ASSERT (mValue<=mMax);
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//             ----  o     ---            ----                               //
//             |   )    |   |    _    |  |      ___    _    ___              //
//             |---  | -+-  |  |/ \  -+- | --- /   ) |/ \  /   )             //
//             |   ) |  |   |  |   |  |  |   \ |---  |   | |---              //
//             |___  |   \ _|_ |   |   \ |___/  \__  |   |  \__              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

BitIntGene::BitIntGene (const GeneticID& id, int mi, int ma, int bits, const StringMap& params, double mut) : AnyIntGene (id, mi, ma, mut) {
	ASSERT (bits>=1 && bits<=32);
	ASSERTWITH (ma-mi+1==(1<<bits), format("Range (%d-%d) must match the number of bits (%d)",
											 mi, ma, bits));
	
	mBitCount = bits;
	for (int i=0; i<bits; i++)
		mBits.add (new BinaryGene (format("b%d", i)));

	Gentainer* dummy = NULL;
	mBits.addPrivateGenes (*dummy, params);
	mBits.selfadjust (false);

	trywith (mGrayCoded = getOrDefault (params, "BitIntGene.graycoding", String(1)).toInt (),
			 (CONSTR) format ("Gene name='%s'", (CONSTR) id));
}

void BitIntGene::copy (const Genstruct& o) {
	AnyIntGene::copy (o);
	shallowCopy (static_cast<const BitIntGene&>(o));
}

void BitIntGene::shallowCopy (const BitIntGene& o) {
	mBitCount	= o.mBitCount;
	mBits.copy (o.mBits);
	mGrayCoded = o.mGrayCoded;
}

int BitIntGene::getvalue () const {
	// Read the bits to a binary vector
	PackArray<int> bits (mBitCount);
	for (int i=0; i<mBitCount; i++)
		bits[i] = static_cast<const BinaryGene&> (*mBits.getGene(format("b%d",i))).getvalue();

	if (mGrayCoded) {
		// Convert from Gray to binary
		for (int i=mBitCount-2; i>=0; i--)
			bits[i] = bits[i+1]^bits[i];
	}

	// Convert from binary to double
	unsigned long int sum = 0;
	for (unsigned long int i=0, m=1; i < (uint) mBitCount; i++, m<<=1)
		sum |= m*bits[i];

	return mMin + sum;
}

void BitIntGene::print (TextOStream& out) const {
	out.printf ("%s=", (CONSTR) id);
	mBits.print (out);
}

void BitIntGene::recombine (const Genstruct& ar, const Genstruct& br) {
	const BitIntGene& a = static_cast<const BitIntGene&> (ar);
	const BitIntGene& b = static_cast<const BitIntGene&> (br);
	
	mBits.recombine (a.mBits, b.mBits);
}

void BitIntGene::check () const {
	AnyIntGene::check ();
	ASSERT (getvalue()>=mMin);
	ASSERT (getvalue()<=mMax);
	mBits.check ();
	ASSERT (mBitCount>=1);
	ASSERT (mBitCount<1000); // Reasonable maximum
}


/*


BitIntGene::BitIntGene (const GeneticID& nam, int mi, int ma, double mut) : Gentainer (nam), IntGene (nam, mi, ma, mut) {
	// Calculate the number of bits
	int range = ma-mi+1;
	mBits=1;
	for (int i=0; i<8 && range>0; i++)
		range = range >> 1;
	
	ASSERTWITH ((1<<mBits)==range, "BitIntGene range has to be power of 2");

	// Add the bit genes
	for (int i=0; i<mBits; i++)
		add (new BinaryGene (format("b%d", i)));

	init ();
}

BitIntGene::BitIntGene (const BitIntGene& o) : IntGene(o), Gentainer(o) {
	mBits = o.mBits;
}

int BitIntGene::getvalue () const {
	int sum=0;
	bool v;
	for (int i=0; i<mBits; i++) {
		const BinaryGene* bit = static_cast<const BinaryGene*>(Gentainer::getGene(format("b%d", i)));
		ASSERT (bit);
		if (bit->getvalue())
			sum += 1<<i;
	}
	return sum+mMin;
}

void BitIntGene::init () {
	Gentainer::init ();
}

void BitIntGene::copy (const Genstruct& o) {
	Gentainer::copy (o);
	IntGene::copy (o);
	const BitIntGene& other = dynamic_cast<const BitIntGene&>(o);
	mBits = other.mBits;
}
	
Genstruct* BitIntGene::replicate () const {
	return static_cast<Genstruct*>(new BitIntGene (*this));
}

bool BitIntGene::pointMutate (const MutationRate& mut_rate) {
	return Gentainer::pointMutate (mut_rate);
}

double BitIntGene::equality (const Genstruct& o) const {
	return Gentainer::equality (o);
}

void BitIntGene::print (TextOStream& out) const {
	out.printf ("%s=%d", (CONSTR) Gentainer::id, getvalue());
}

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//              ---                      ----                                //
//               |    _    |   ___      |      ___    _    ___               //
//               |  |/ \  -+- /   ) |/\ | --- /   ) |/ \  /   )              //
//               |  |   |  |  |---  |   |   \ |---  |   | |---               //
//              _|_ |   |   \  \__  |   |___/  \__  |   |  \__               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool InterGene::execute (const GeneticMsg& msg) const {
	return msg.mrHost.execute (GeneticMsg (targetGene, msg.mrHost));
}

void InterGene::copy (const Genstruct& o) {
	copyGenstr (o);
	shallowCopy (static_cast<const InterGene&>(o));
}	

void InterGene::print (TextOStream& out) const {
	out.printf ("%s->%s", (CONSTR) id, (CONSTR) targetGene);
}
