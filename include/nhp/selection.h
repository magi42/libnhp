#ifndef __SELECTION_H__
#define __SELECTION_H__

#include <magic/mmatrix.h>
#include <magic/mrefarray.h>
#include "nhp/individual.h"

// Externals
class SimplePopulation;
class Gentainer;

class SelectionSituation;


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//    ----       |                 o            |   |               o       //
//   (      ___  |  ___   ___   |           _   |\ /|  ___   |              //
//    ---  /   ) | /   ) |   \ -+- |  __  |/ \  | V |  ___| -+- |/\ | \ /   //
//       ) |---  | |---  |      |  | /  \ |   | | | | (   |  |  |   |  X    //
//   ___/   \__  |  \__   \__/   \ | \__/ |   | |   |  \__|   \ |   | / \   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** Matrix for calculating the selection probabilities of specific
 *  parent pairs.
 *
 *  Yes, the selection mechanism acts on parent PAIRS; it doesn't just
 *  select one parent and then another parent and then mate them. This
 *  mechanism allows individuals themselves to affect their
 *  selection. It is also easy to incorporate various other factors
 *  such as distance in some environment-space as a selection
 *  probability factor between two individuals.
 **/
class SelectionMatrix : public Object {
  public:

	/** Standard constructor.
	 *
	 **/
	explicit	SelectionMatrix		(const SelectionSituation& situation);
				~SelectionMatrix	() {}
	
	/** Selects a random pair of individuals from the
	 *  population. These individuals will be mated.
	 *
	 *  @param a Ordered (!) population index number of the first, selecting parent.
	 *  @param b Ordered (!) population index number of the first, selecting parent.
	 **/
	void	selectRandomPair		(int& a, int& b) const;
	
	/** Implementation for @ref Object. */
	TextOStream&	operator>>		(TextOStream& out) const {return mSelection >> out;}

  protected:

	/** Calculates the selection matrix.
	 **/
	void	calculateMatrix			(const SelectionSituation& situation);

  private:
	SelectionMatrix();

	/** Resulting selection matrix */
	Matrix	mSelection;
	
	friend class Selector;
};



/////////////////////////////////////////////////////////////////////////////////////////
//  ----       |                 o             ---- o                     o            //
// (      ___  |  ___   ___   |           _   (        |         ___   |           _   //
//  ---  /   ) | /   ) |   \ -+- |  __  |/ \   ---  | -+- |   |  ___| -+- |  __  |/ \  //
//     ) |---  | |---  |      |  | /  \ |   |     ) |  |  |   | (   |  |  | /  \ |   | //
// ___/   \__  |  \__   \__/   \ | \__/ |   | ___/  |   \  \__!  \__|   \ | \__/ |   | //
/////////////////////////////////////////////////////////////////////////////////////////

/** A temporary situation where pairs of are expected to be select for
 *  mating.
 *
 *  Handles the selection calculations for a population's individuals;
 *  generates a @ref SelectionMatrix.
 *
 *  This class is primarily used by @ref EAStrategy and similar classes.
 **/
class SelectionSituation : public Object {
  public:

	/** Standard constructor.
	 *
	 *  @param pop The population.
	 *
	 *  @param ordpop Individuals of the population ordered according
	 *  to their fitnesses.
	 **/
							SelectionSituation		(const SimplePopulation& pop);

	const SimplePopulation&	population				() const {return mrPop;}

	/** Returns the i:th fittest @ref Individual from the
	 *  population. Fittest individual has index 0! Const version.
	 **/
	const Individual&		getOrdered				(int i) const {return mOrdPop[i];}
	
  private:

	/** Actual population */
	const SimplePopulation& mrPop;

	/** Ordered SimplePopulation */
	RefArray<Individual> mOrdPop;
};



//////////////////////////////////////////////////////////////////////////////
//      ----       |                 o            ----                      //
//     (      ___  |  ___   ___   |           _   |   )            ____     //
//      ---  /   ) | /   ) |   \ -+- |  __  |/ \  |---  |/\ |/|/| (         //
//         ) |---  | |---  |      |  | /  \ |   | |     |   | | |  \__      //
//     ___/   \__  |  \__   \__/   \ | \__/ |   | |     |   | | | ____)     //
//////////////////////////////////////////////////////////////////////////////

/** A container object that contains selection parameters for
 *  different selection methods.
 *
 *  There is typically exactly one of these objects per population,
 *  owned by @ref Population as global selection parameters.
 **/
class SelectionPrms : public Object {
  public:

	explicit				SelectionPrms	();
	explicit				SelectionPrms	(const SelectionPrms& o) {copy(o);}
							~SelectionPrms	() {}

	/** Sets the selection propability of one selection method to 1.0 and others to 0.0.
	 **/
	void					useOnlyMethod	(int m);

	/** Set the self-adaptation of different selection method parameters
	 **/
	void					adaptParams		(bool mu=true, bool etaPlus=true,
											 bool q=true);

	void					adaptWeights	(bool v=true) {mAdaptiveWeights=v;}

	/** Sets the mu-parameter for (mu,+lambda)-Selection.
	 **/
	void					setMu		(int m);
	
	/** Sets the mu-parameter for (mu,+lambda)-Selection as a
	 *  partial value [0,1] of the population size.
	 **/
	void					setMuPart	(double part);

	/** Sets the etaPlus parameter for proportional selection.
	 **/
	void					setEtaPlus		(double etaplus);

	/** Sets the q parameter for tournament selection.
	 **/
	void					setQ			(int q);

	/** Returns mu for specified population size (required since
	 *  micro may have been given as a percentage).
	**/
	int						muFor			(int populsize) const;
	
	void					copy			(const SelectionPrms& o);

  protected:
	/** mu-parameter for (mu,+lambda)-Selection as a portion of the
	 *  population size. Domain is (0,1).
	 **/
	double					mMuPart;

	/** mu-parameter for (mu,+lambda)-Selection. */
	int						mMu;

	/** Mode flag: shall we use self-adaptive mu parameter? */
	bool					mAdaptiveMu;

	/** Selection parameter for proportional selection. */
	double					mEtaPlus;

	/** Mode flag: shall we use self-adaptive etaPlus parameter? */
	bool					mAdaptiveEtaPlus;

	/** Selection parameter for tournament selection. */
	int						mQ;

	/** Mode flag: shall we use self-adaptive q parameter? */
	bool					mAdaptiveQ;

	/** Mode flag: Should we weigh between different selection
	 *  schemes (true) or should we select one method with the
	 *  propability given by it's weight (false). Default=false (it's
	 *  faster).
	 **/
	bool					mWeightedSelection;

	/** Weights OR probabilities of different selection methods.
	 **/
	Vector					mSelMethodW;

	/** Mode flag: are the selection method weights or probabilities
	 *  self-adaptive?
	 **/
	bool					mAdaptiveWeights;

	friend class Selector;
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                   ----       |                                           //
//                  (      ___  |  ___   ___   |                            //
//                   ---  /   ) | /   ) |   \ -+-  __  |/\                  //
//                      ) |---  | |---  |      |  /  \ |                    //
//                  ___/   \__  |  \__   \__/   \ \__/ |                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** Selection handler of an @ref Individual.
 *
 *  Each individual has one instance of this class to handle all it's
 *  selection tasks.
 *
 *  This class provides some additional functionality for @ref
 *  SelectionPrms. It is kept separate because we don't want to give
 *  these methods if not necessary.
 **/
class Selector : public SelectionPrms {
  public:
	explicit				Selector		() : SelectionPrms () {mScore=0;}
	explicit				Selector		(const SelectionPrms& orig);
							~Selector		() {}
	
	/** Inserts genes that control selection into the given genome, if
	 *  the selection parameters are set to be self-adaptive.
	 *
	 *  @param params Dynamic parameters in a @ref String @ref
	 *  Map. See various methods of @ref SelectionPrms and @ref
	 *  SelectionMatrix for further descriptions of the parameters.
	 *
	 *  @param params["mu"] mu-parameter for (mu,+lambda)-Selection.
	 *  @param params["q"] q-parameter for tournament selection.
	 *  @param params["eta+"] eta+-parameter for proportional selection.
	 *  @param params["adaptMu"] Should mu be self-adapted? [Default:0]
	 *  @param params["adaptQ"] Should q be self-adapted? [Default:0]
	 *  @param params["adaptEta+"] Should eta+ be self-adapted? [Default:0]
	 **/
	static void				addGenesTo		(Gentainer& g, const StringMap& params);

	/** Reads selection parameters from the given genome.
	 **/
	void					read			(const Genome& g);

	/** Returns the selection affinity of individual i to individual
	 *  j, combined over all different selection methods.
	 **/
	double					select			(const SelectionSituation& situation,
											 int self_i, int other_j) const;

	/** Tells how many times the individual has got some. (Hmm, shouldn't
	 *  this kind of information be private??).
	**/
	int						score			() const {return mScore;}
	
	/** Makes the individual virgin again (resets it's score). Wow. **/
	void					virginize		() {mScore=0;}

	enum smconsts {MULAMBDASELECTION=0, LINEARRANKING, PROPORTIONAL, TOURNAMENT,
				   number_of_methods};
	
  protected:
	/** Returns the selection probability for individuals i and j.
	 *
	 *  @param smno Selection method number (see the enum 'smconsts').
	 *  @param i The individual making the selection choise.
	 *  @param j The target of selection.
	 *
	 *  @return Selection probability (not affinity) of Individual j by i.
	 **/
	virtual double	selectWithMethod	(const SelectionSituation& situation, int smno, int i, int j) const;

	/** Evolution Strategies (mu,+lambda)-Selection, where lambda is
	 *  the number fittest individuals in the population that can be
	 *  potential parents, and lambda is the number of offspring.
	 **/
	double	muLambdaSelection			(const SelectionSituation& situation, int i, int j) const;

	/** Selection method where the individuals are ranked according to
	 *  their fitness, and the selection probability is a linear
	 *  function of this rank.
	 **/
	double	linearRanking				(const SelectionSituation& situation, int i, int j) const;

	/** Proportional selection.
	 **/
	double	proportionalSelection		(const SelectionSituation& situation, int i, int j) const;

	/** Tournament selection.
	 **/
	double	tournamentSelection			(const SelectionSituation& situation, int i, int j) const;

  protected:
	/** Some selection methods, like the method the canonic GA, select
	 *  some of the best individuals with 100% propability. To achieve
	 *  this, we have to count how many times an individual has got
	 *  some.
	 **/
	int						mScore;

	void					operator=		(const Selector& other) {}
	void					operator=		(const SelectionPrms& other) {}
};

/*
class SelectionMethod : public Object {
  public:
						SelectionMethod	(const String& name);

	const Array<SelectionMethod>& methods	() const {return sMethods;}
	virtual	double			select			(int i, int j);
	
  private:
	static Array<SelectionMethod> sMethods;
};
*/

#endif
