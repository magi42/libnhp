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

#include <ctype.h>
#include <magic/mmath.h>
#include <magic/mstream.h>
#include <magic/mclass.h>

#include "nhp/genetics.h"
#include "nhp/individual.h"
#include "nhp/gaenvrnmt.h"

impl_abstract (EAEnvironment, {Object});



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//    ----   _   -----             o                                        //
//   |      / \  |       _                      _          ___    _    |    //
//   | --- /   \ |---  |/ \  |   | | |/\  __  |/ \  |/|/| /   ) |/ \  -+-   //
//   |   \ |---| |     |   |  \ /  | |   /  \ |   | | | | |---  |   |  |    //
//   |___/ |   | |____ |   |   V   | |   \__/ |   | | | |  \__  |   |   \   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

EAEnvironment::EAEnvironment ()
{
	mNEvals			= 1;
	mNoise			= 0.0;
	bestfitn		= 1E+30;
	mBestFitness	= 1E+30;
	mTotEvals		= 0;
	mpBest			= NULL;
	mCycles			= 0;
}

/*******************************************************************************
* Evaluates fitness of an individual.
*
* If artificial noise is defined for the environment, applies it to the
* fitness.
*
* Keeps record of the best measured fitness and stores the individual
* having the fitness until a better one is found.
*
* @return Measured fitness.
*******************************************************************************/
double EAEnvironment::evaluate (const Individual& ind)
{
	// Evaluate the fitness
	double fitness = evaluateg (ind);

	// Record the best _objective_ fitness. Note that this is done
	// before adding the artificial noise, so this is really the true
	// fitness.
	if (fitness < bestfitn)
		bestfitn = fitness;

	// Add some artificial noise. The if statement is here because we don't
	// want to compate to 0.0.
	if (mNoise > 0.0001)
		fitness += gaussrnd (mNoise);

	// Record the best _subjective_ fitness
	if (fitness < mBestFitness) {
		mBestFitness = fitness;
		mpBest       = const_cast <Individual*> (&ind);
	}

	mTotEvals++;

	return fitness;
}

void EAEnvironment::check () const
{
	ASSERT (mNEvals>=1 && mNEvals<100000);
	ASSERT (mTotEvals>=0 && mTotEvals<1000000);
	ASSERT (mCycles>=0 && mCycles<100000);
	ASSERT (mNoise>=0.0 && mNoise<10000.0);
}
