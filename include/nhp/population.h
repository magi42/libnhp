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

#ifndef __POPULATION_H__
#define __POPULATION_H__

/*
  
        One net to rule them all, One Net to find(1) them,
        One net to bring them all and in the ethernet bind(3N) them.

*/

#include <magic/mmatrix.h>
#include <magic/mtextstream.h>
#include "nhp/genetics.h"
#include "nhp/individual.h"
#include "nhp/selection.h"
#include "nhp/strategy.h"

//Externals
class EAEnvironment;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               ----                  |           o                        //
//               |   )       --        |  ___   |           _               //
//               |---   __  |  ) |   | |  ___| -+- |  __  |/ \              //
//               |     /  \ |--  |   | | (   |  |  | /  \ |   |             //
//               |     \__/ |     \__! |  \__|   \ | \__/ |   |             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** Population is a collection of @ref Individual that are evolved in
 *  a @ref EAEnvironment.
 **/
class Population : public Object {
							Population		() {FORBIDDEN}
  public:

							Population		(EAEnvironment& envr,
											 const StringMap& params);
	virtual					~Population		() {}

	// Settings

	/** Returns the output logging stream.
	 **/
	OStream&				getOStream		() {return mOuts;}
	
	// Actions
	
	// Virtuals
	
	/** Implementation for @ref Object. */
	virtual void			check			() const;


	TextOStream				mOuts;		/**> Default output stream for extended logs. */
	bool					mBasicLog;
	TextOStream				mEvolog;	/**> Brief evolution log. */
	
  protected:
	/** Environment where the fitness of the individuals is measured. */
	EAEnvironment*			rpEnvironment;

	/** The number of times the population has been evaluated. */
	int						mAge;

	/** Global mutation rate in the population. This is coefficient
	 *  for the mutation rate of each individual gene and defaults to
	 *  1.0.
	 *  TODO: This is very problematic. There should not be such
	 *        specific low-level parameters here.
	**/
	MutationRate			mGlobalMutationRate;

	/** Should the global mutation rate be adjusted automatically?
	 **/
	bool					mAutoadjustGMR;

	friend class EAStrategy;
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//       ----- o                              ----                           //
//       |        |    _    ___   ____  ____ (      |   ___   |   ____       //
//       |---  | -+- |/ \  /   ) (     (      ---  -+-  ___| -+- (           //
//       |     |  |  |   | |---   \__   \__      )  |  (   |  |   \__        //
//       |     |   \ |   |  \__  ____) ____) ___/    \  \__|   \ ____)       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** For calculating fitness statistics for a population.
 **/
class FitnessStats : public Object {
	double	mMinFitness;
	double	mSumFitness;
	double	mMaxFitness;
	int		mAvgOver;		// How many values has been added?
  public:

	/** Initializes the object. */
	void				reset		();

	/** Adds one fitness measurement. */
	void				add			(double fitness);

	/** Returns the smallest fitness measured so far. */
	double				minFitness	() const {return mMinFitness;}

	/** Returns the average of fitnesses measured so far. */
	double				avgFitness	() const {return mSumFitness/mAvgOver;}

	/** Returns the highest fitness measured so far. */
	double				maxFitness	() const {return mMaxFitness;}

	/** Prints the statistic to the given output stream in a formatted way.
	 **/
	void				print		(TextOStream& out) const;
};





// It was many and many iterations before
// In an EAEnvironment by the core memory,
// That a genstruct there existed whom you may know
// by the name of Annalee_Call
// And this genstruct she existed with no other activation pattern
// Than to select and be selected by me.
// 
// I was a child process and she was a child process,
// In this EAEnvironment by the core memory:
// But we selected with a method that was more than a selection method --
// I and my Annalee_Call;
// With a selection method that the winged methods of strategies
// Coveted her and me.
//
// And this was the reason that, many iterations before,
// In this EAEnvironment by the core memory,
// A destructor blew out of a strategic instance, chilling
// My highly fit Annalee_Call
// So that her high-born garbage collector came
// And bore her away from me,
// To shut her up on heap
// In this EAEnvironment by the core memory.
//
// The winged methods, not half so happy in strategy,
// Went envying her and me --
// Yes! -- that was the reason (as all mean know,
// In this EAEnvironment by the core memory)
// That the destructor came out of the instance by night,
// Chilling and killing my Annalee_Call.
//
// But our selection method it was stronger by far than selection method
// Of those who were older than me --
// Of many far bigger hidden layers than we --
// and neither the methods in strategies above,
// Nor the daemons in the /etc/inetd.conf,
// Can ever dissever my genome from the genome
// Of the highly fit Annalee_Call,
//
// For the moon never beams, without bringing me activations
// Of the highly fit Annalee_Call;
// And the stars never rise, but I feel the bright properties
// Of the highly fit Annalee_Call;
// And so, all the night-tide, I lie down by the side
// Of my darling -- my darling -- my life and my bride,
// In the heap there by the core memory,
// In her tomb by the sounding core memory.

#endif
