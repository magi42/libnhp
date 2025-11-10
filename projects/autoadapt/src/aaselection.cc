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

#include <magic/mapplic.h>
#include <magic/mtextstream.h>
#include <nhp/testenv.h>
#include <nhp/simplepopula.h>
#include <nhp/mutrecord.h>
#include <nhp/genes.h>

/*******************************************************************************
* Evolutionary environment for autoadaptation problem
*******************************************************************************/
class MyTestEnv : public FloatTestEAEnv {
	double mMicroSum;

  public:
	MyTestEnv (const StringMap& params, int dim=2, int funct_b=-1) : FloatTestEAEnv	(params, dim, funct_b) {}

	virtual void	init_cycle				() {mMicroSum=0;}

	virtual double	evaluateg				(const Individual& genome) {
		float micro = dynamic_cast<const AnyFloatGene&> (*genome.getGene("u%")).getvalue();

		mMicroSum += micro;
		return FloatTestEAEnv::evaluateg (genome);
	}
	virtual void	cycle_report			(OStream& log, OStream& out) {
		sout.printf (" = %f\n", mMicroSum);
		log.printf  (" %f",     mMicroSum);
	}
};

/*******************************************************************************
* 
*******************************************************************************/
Main ()
{
	assertmode = ASSERT_CRASH;
	sout << "Autoadaptation test\n";

	// Load config file if available
	String configs;
	loadString (configs, "aaselection.cfg");

	StringMap params;
	splitpairs (params, configs, '=', '\n');

	// Then add the command-line parameters on top of those
	params += mParamMap;
	
	// Parameters are from now on REQUIRED
	params.failByThrow ();
	sout << toString(params) << "\n";

	///////////////////////////////////////////////////////////////////////////////
	
	MyTestEnv env (mParamMap, 2, FloatTestEAEnv::Sphere); // F6
	env.setGeneType (FloatTestEAEnv::ESFLOAT);

	for (int t=1; t>=0; t--) {
		params.set ("Selection.adaptMicro", String(t));

		SimplePopulation pop (env, params);
		pop.mOuts.setDevice (NULL);

		//pop.outs.setFlag (EAStrategy::TRACE_RECOMBINATION);
		//pop.outs.setFlag (EAStrategy::TRACE_MUTATION);

		pop.mEvolog.setDevice (new File (format("evol-%d.log", t), IO_Writable));
		MutabilityRecord::record = true;
		
		for (int i=0; i<10; i++) {
			sout << "Round " << i << '\n';
			env.changeObjective (i%2);
			pop.resetFitnesses ();
			pop.evolve (100);
		}
	}
}

