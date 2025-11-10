#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include "nhp/genetics.h"
#include <magic/mmath.h>
#include <magic/mmap.h>

//Externals
class EAEnvironment;
class SimplePopulation;
class SelectionMatrix;
class Selector;
class SelectionPrms;



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//               ---           | o       o     |             |               //
//                |    _       |               |        ___  |               //
//                |  |/ \   ---| | |   | |  ---| |   |  ___| |               //
//                |  |   | (   | |  \ /  | (   | |   | (   | |               //
//               _|_ |   |  ---| |   V   |  ---|  \__!  \__| |               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** The carrier of a genome; can have phenotype.
 *
 *  Another central feature of an @ref Individual is that it's fitness
 *  in an @ref EAEnvironment can be measured.
 **/
class Individual : public Comparable {
	decl_dynamic (Individual);
  public:

							Individual		();
	explicit				Individual		(const Individual& prototype);

	/** Construction from genotypic template.
	 **/
	explicit				Individual		(const Genome& prototype);

	virtual					~Individual		();
	
	/** Resets the individual to birth conditions; removes any
	 *  phenotypic features.
	 *
	 *  @param init Giving a TRUE value causes the genome of
	 *  the individual to be initialized (randomized).
	**/
	void					incarnate		(bool init);

	///////////////////////////////////////////////////////////////////////////////
	// Passthroughs to feature map
	
	/** Add or change a phenotypic feature. */
	void					set				(const String& key, Object* val) {
		mFeatures.set (key, val);
	}

	/** Remove a phenotypic feature. */
	void					remove			(const String& key) {mFeatures.remove(key);}
	
	/** Returns a phenotypic feature. */
	Object&					operator[]		(const String& key) {return mFeatures[key];}

	/** Returns a const phenotypic feature. */
	const Object&			operator[]		(const String& key) const {return mFeatures[key];}

	Object*					getFeature		(const String& key) {return mFeatures.getvp(key);}
	const Object*			getFeature		(const String& key) const {return mFeatures.getp(key);}

	/** Returns a reference to the phenotypic features of the
	 *  Individual.
	 **/
	const Map<String,Object>&	features	() const {return mFeatures;}
	
	///////////////////////////////////////////////////////////////////////////////
	// Evaluation of the fitness
	
	double					evaluate		(EAEnvironment& envr, bool force=false);

	// Lets the individual to give it's preference for the given individual
	//double					select			(const RefArray<Individual>& opop, int j) const;

	// We want to be comparable by our fitness to allow sorting
	int						compare			(const Comparable& other) const;

	// Returns the evaluated fitness of the phenotype
	double					getfitness		() const {return fitness;}

	// Returns the number of evaluations the calculated fitness has
	// been averaged over
	int						averaged_over	() const {return avg_over;}

	// Joins the fitness with a (practically) identical other phenotype
	void					joinfitness		(const Individual& other);

	// Reset fitness and it's averaging
	void					resetFitness	() {fitness=0; avg_over=0;}
	
	// Feature shortcuts

	/** Returns the age of the individual; how many generations it has
	 *  been an elite.
	 **/
	int						getage			() const {return age;}

	/** Returns the selection handler of the Individual.
	 **/
	const Selector&			selector		() const {return *mpSelector;}

	/** Sets the selector of the Individual using the given
	 *  template. This method is typically called by a @ref Population
	 *  or it's @ref EAStrategy, and the template consists of global
	 *  selection parameters.
	 **/
	void					setSelector		(const SelectionPrms& templ);
	
	/** Tells how much the individual likes another.
	 *
	 *  @param selmat Selection matrix.
	 *  
	 *  @return Selection affinity.
	 **/
	//double					select			(const SelectionMatrix& selmat,
	//										 int self_i, int other_j) const;

	// Passthroughs to Genome

	/** Passthrough to the @ref Genome of the Individual. */
	bool				execute				(const GeneticMsg& msg) const {return genome.execute(msg);}
	/** Passthrough to the @ref Genome of the Individual. */
	void				init				() {genome.init();}
	/** Passthrough to the @ref Genome of the Individual. */
	bool				pointMutate			(const MutationRate& k);
	/** Passthrough to the @ref Genome of the Individual. */
	const Genstruct*	getGene				(const GeneticID& n) const {return genome.getGene(n);}
	/** Passthrough to the @ref Genome of the Individual. */
	static void			addGenesTo			(Genome& g, const StringMap& params);
	/** Passthrough to the @ref Genome of the Individual. */
	void				addking				() {genome.addking();}
	/** Passthrough to the @ref Genome of the Individual. */
	int					getkings			() const {return genome.getkings();}
	/** Passthrough to the @ref Genome of the Individual. */
	void				recombine			(const Individual& a, const Individual& b);
	/** Passthrough to the @ref Genome of the Individual. */
	double				equality			(const Individual& other) const {
		return genome.equality (other.genome);
	}
	
	/** Brief printout. */
	void					print			(TextOStream& out) const;

	/** Implementation for @ref Object. Verbose printout. */
	virtual DataOStream&	operator>>		(DataOStream& out) const;

	/** Implementation for @ref Object. */
	virtual void			check			() const;

  private:
	/** Dummy. */
	int						operator==		(const Comparable& other) const {FORBIDDEN}

	/** Adds 1 to age. */
	void					grow_older		();

	/** The phenotypical features of the specimen.
	 **/
	Map<String,Object>		mFeatures;

	/** The average cached fitness of the equal phenotypes of the genome. */
	mutable double			fitness;

	/** From how many evaluations this fitness has been averaged. */
	mutable int				avg_over;

	/** Even artificial lifeforms get older. Isn't it sad? */
	int						age;

	/** Selection handler. */
	Selector*				mpSelector;

	/** The genome of the individual. */
	Genome					genome;

	Individual& operator= (const Individual& other) {FORBIDDEN; return *this;} // Prevent copying
};

#endif
