#include <magic/mapplic.h>
#include <magic/mtextstream.h>
#include <nhp/testenv.h>
#include <nhp/metapopulation.h>
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

