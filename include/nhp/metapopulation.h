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

#ifndef __METAPOPULATION_H__
#define __METAPOPULATION_H__

#include <population.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    |   |                                      |           o               //
//    |\ /|  ___   |   ___   --        --        |  ___   |           _      //
//    | V | /   ) -+-  ___| |  )  __  |  ) |   | |  ___| -+- |  __  |/ \     //
//    | | | |---   |  (   | |--  /  \ |--  |   | | (   |  |  | /  \ |   |    //
//    |   |  \__    \  \__| |    \__/ |     \__! |  \__|   \ | \__/ |   |    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
 * Base class for metapopulation models.
 *****************************************************************************/
class Metapopulation : public Population {
  public:
				Metapopulation (EAEnvironment& envr,
								const StringMap& params)
						: Population (envr, params) {
				}

	virtual		~Metapopulation (){}

  protected:
};

/******************************************************************************
 * Metapopulation with uniform probability for migration from any
 * patch to another patch.
 *****************************************************************************/
class NonspatialMetapopulation : public Metapopulation {
  protected:
};

#endif
