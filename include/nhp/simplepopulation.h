#include <population.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  ---- o            |       ----                  |           o            //
// (              --  |  ___  |   )       --        |  ___   |           _   //
//  ---  | |/|/| |  ) | /   ) |---   __  |  ) |   | |  ___| -+- |  __  |/ \  //
//     ) | | | | |--  | |---  |     /  \ |--  |   | | (   |  |  | /  \ |   | //
// ___/  | | | | |    |  \__  |     \__/ |     \__! |  \__|   \ | \__/ |   | //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class SimplePopulation : public Population {
	// The current set of individuals in the population
	Array<Individual>*		population;

	// Some statistics of the current population
	FitnessStats			mFitnessStats;

	// Number of elites, individuals who will survive intact to the
	// next generation
	int						mElites;
	bool					mUseGlobalElites;

	SelectionPrms			mSelectionParams;
	
	// Global micro portion The number of potential parents. The
	// semantics are dependent on the evolutionary strategy used
	bool					mUseGlobalMicro;

	// Global q parameter for tournament selection
	bool					mUseGlobalQ;

	// Global etaPlus parameter for linear ranking selection
	bool					mUseGlobalEtaPlus;
	
	// Global mutation rate in the population. This is coefficient for
	// the mutation rate of each individual gene and defaults to 1.0
	MutationRate			mGlobalMutationRate;

	// Should the global mutation rate be adjusted automatically?
	bool					autoadj_gmr;

  public:

	// Minimum allowed similarity between genomes
	double					minsimilarity;
	
	// Constructs a binary population of the given size
							SimplePopulation		(EAEnvironment& envr,
											 const StringMap& params);

							~SimplePopulation		();
	
	double					evolve			(int generations, const char* savefile=NULL,
										 double trg_fitn=-1);

	// Returns an individual in the population by it's index number
	const Individual&		operator[]		(int i) const {return (*population)[i];}

	// Same non-const
	Individual&				operator[]		(int i) {return (*population)[i];}
	
	// Returns the number of individuals in the population
	int						size			() const {return population->size;}

	// Dumps the population to given output
	void					print			(DataOStream& out) const;

	// Writes an one-line generation report to the given stream
	void					report			(OStream& log) const;

	// Returns the current strategy
	const EAStrategy&		getstrategy		() const {return *strategy;}

	// Returns the number of times the population has been evaluated
	// (i.e., the generation)
	int						getAge			() const {return mAge;}

	// Set the global mutation rate
	MutationRate&			mutRate			() {return mGlobalMutationRate;}

	SelectionPrms&			selParams		() {return mSelectionParams;}
	const SelectionPrms&	selParams		() const {return mSelectionParams;}

	// Resets the fitness averaging of individuals. Useful when
	// changing the objective function (old fitness values would be
	// invalid).
	void					resetFitnesses	();

	void					check			() const;
	
  private:

	// Evaluates the population in an environment
	void					evaluate		(EAEnvironment& envr, OStream& out);

 	// Add population-dependent features to a genome. I suppose there
	// might be some use for this. Maybe.
	void					addFeaturesTo	(Genome& genome) const;

	friend EAStrategy;
	/*
	friend CommonStrategy;
	friend CommaStrategy;
	friend PlusStrategy;
	*/
	friend SelectionMatrix;
};
