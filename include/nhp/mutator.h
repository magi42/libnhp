#ifndef __MUTATOR_H__
#define __MUTATOR_H__

#include <magic/mmath.h>

/** Interface for @ref FloatGene mutation operators.
 **/
class FloatMutator {
  public:
	virtual ~FloatMutator () {}
	
	/** Mutates the given current value with the mutation
	 *  function. All inheritors should implement this.
	 *
	 *  @param x The current value of the gene.
	 *  @param min Minimum value of the gene.
	 *  @param max Maximum value of the gene.
	 *
	 *  @param variance Distribution parameter for the mutation. It
	 *  can be variance or standard deviation for normally distributed
	 *  mutations, but it's semantics are really dependent on the
	 *  actual mutation method.
	 **/
	virtual double	mutate	(double x, double min, double max, double variance) const{FORBIDDEN}
};

/** The standard normally distributed @ref FloatGene mutator. It
 *  ignores the variance given in the parameter. This is useful for
 *  protecting specific genes from the self-adaptation of mutation
 *  rates.
 **/
class ConstantMutator : public FloatMutator {
	double mVariance;
  public:
	/** Standard constructor.
	 *
	 *  @param var Variance. This value is IGNORED! Instead, the
	 *  variance given as parameter to the mutate () method is used.
	 **/
					ConstantMutator	(double var=0.22) : mVariance(var){}

	/** As in the superclass, except that the variance parameter is
	 *  not used here; instead the value given in the constructor is used.
	 **/
	virtual double	mutate	(double p, double min, double max, double variance) const {
		//TRACE2 ("p=%f, v=%f", p, variance);
		double pp=p;
		if (max>min) { // Can be 0 -> gene is immutable
			// Default mutation
			double delta;
			do {
				delta = gaussrnd (mVariance);
			} while (pp+delta<min || pp+delta>max);
			pp += delta;
		}
		return pp;
	}
};

/** @ref FloatGene mutator for mutating mutation rates.
 **/
class MutationMutator : public FloatMutator {
	double mPhi;
  public:
	/** Standard constructor.
	 *
	 *  @param phi Distribution parameter for the method.
	 **/
					MutationMutator	(double phi=0.22) : mPhi(phi){}
	
	/** As in the superclass, except that the variance parameter is
	 *  not used here; instead the phi value given in the constructor is used.
	 *
	 *  @param variance NOT USED
	 **/
	virtual double	mutate	(double p, double min, double max, double variance) const {
		//TRACE2 ("p=%f, v=%f", p, variance);
		double pp = 1/(1+(1-p)/p*exp(-gaussrnd(mPhi)));
		if (pp<min)
			pp=min;
		return pp;
	}
};

/** @ref FloatGene mutator from Thomas Bäck 1997, SNAC Book, page 37 (10).
 **/
class MutationSimpleMutator : public FloatMutator {
	double mTau;
  public:
	/** The standard constructor.
	 *
	 *  @param n Number of genes in the system. In pure Evolution
	 *  Strategies with only real-valued genes this would mean the
	 *  length of the real vector, i.e. the number of real genes in
	 *  the genome.
	 **/
					MutationSimpleMutator	(int n=100) : mTau(1/sqrt(double(n))){}
	
	/** As in the superclass, except that the variance parameter is
	 *  not used here; instead a parameter calculated from the value
	 *  given in the constructor is used.
	 **/
	virtual double	mutate	(double p, double min, double max, double variance) const {
		//TRACE2 ("p=%f, v=%f", p, variance);
		double pp = p*exp(mTau*gaussrnd(1));
		if (pp<min)
			pp=min;
		return pp;
	}
};


/** Normally distributed @ref FloatGene mutator for circular domains.
 **/
class CircularMutator : public FloatMutator {
	double mPhi;
  public:
	/** As in the superclass, except that the variance parameter is
	 *  not used here; instead the value given in the constructor is used.
	 **/
	virtual double	mutate	(double p, double min, double max, double variance) const {
		return 1/(1+(1-p)/p*exp(-mPhi*gaussrnd(1)));
	}
};

#endif
