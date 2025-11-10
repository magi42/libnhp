#ifndef __STRATEGY_H__
#define __STRATEGY_H__

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           -----   _    ----                                               //
//           |      / \  (      |       ___   |   ___                        //
//           |---  /   \  ---  -+- |/\  ___| -+- /   )  ___  \   |           //
//           |     |---|     )  |  |   (   |  |  |---  (   \  \  |           //
//           |____ |   | ___/    \ |    \__|   \  \__   ---/   \_/           //
//                                                      __/   \_/            //
///////////////////////////////////////////////////////////////////////////////

// In olden times strategy was listed among the Ten Abilities and
// Seven Arts as a beneficial practice. It was certainly an art but as
// a beneficial practice it was not limited to sword-fencing. The true
// value of sword-fencing cannot be seen within the confines of
// sword-fencing technique.

/** Search strategy for an evolutionary algorithm. You should not
 *  confuse this with Evolution Strategies, which is one kind of
 *  evolution strategy for real-valued genes.
 *
 *  The strategy is distinguished from the population, as alternative
 *  strategies could be used for a population.
 **/
class EAStrategy {
  public:

	/** Attaches the strategy to the given population. The class
	 *  currently supports only linear @ref SimplePopulation.
	 **/
					EAStrategy		(SimplePopulation& popula);
	virtual			~EAStrategy		() {delete mpNextGen;}

	/** Evolves a population to adapt to an environment
	 *
	 *  @param envr Environment where the population evolves in.
	 *
	 *  @param out Output stream.
	 *
	 *  @param log Logging stream for brief evolution logs.
	 **/
	void			evolve			(EAEnvironment& envr, TextOStream& out, TextOStream& log);
	
	/** Adds strategy-dependent features to the given genome.
	 **/
	virtual void	addFeaturesTo	(Genome& genome) const;

	/** Prints some strategic information to the given stream.
	 **/
	void			print			(TextOStream& out);

	/** Forms the next generation, usually by recombining parent
	 *  individuals as offspring.
	 *
	 *  @param pop_order The population of individuals in an array
	 *  that is ordered according to their fitnesses.
	 *
	 *  @param selmat Selection matrix that contains the pairwise
	 *  selection probabilities for the potential pairs of
	 *  parents. See @ref SelectionMatrix for more information about
	 *  selection.
	 **/
	void			recombine		(const SelectionSituation& situation,
									 const SelectionMatrix& selmat);

	/** @ref OStream flags for outputting some trace information.
	 **/
	enum traceflags {TRACE_RECOMBINATION=10, TRACE_MUTATION};

	/** Implementation for @ref Object */
	void			check			() const;
	
  protected:
	/** The population that is being evolved with this strategy. */
	SimplePopulation&	mrPopula;

	/** An array holding the individuals of "next generation". */
	Array<Individual>*	mpNextGen;

	/** Mode flag dictating whether to allow self-breeding or not. */
	bool allow_same_parents;
};


//
// Timing in strategy
//
// There is timing in everything. Timing in strategy cannot be
// mastered without a great deal of practice.
//
// Timing is important in dancing and pipe or string music, for they
// are in rhythm only if timing is good. Timing and rhythm are also
// involved in the military arts, shooting bows and guns, and riding
// horses. In all skills and ablilities there is timing.
//   -- Miyamoto Musashi
//

#endif
