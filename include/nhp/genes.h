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

#ifndef __GENES_H__
#define __GENES_H__

#include <magic/mobject.h>
#include <magic/mdatastream.h>
#include "nhp/genetics.h"

// Is it a bool? Is it a double? No! It's the SuperGene!

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            ----                                           //
//                           |      ___    _    ___                          //
//                           | --- /   ) |/ \  /   )                         //
//                           |   \ |---  |   | |---                          //
//                           |___/  \__  |   |  \__                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** The abstract concept of a genetic structure encoding some specific
 * phenotypic feature; the atomic unit of a @ref Genstruct.
 *
 * Actually it's not so atomic at physical level, since it can contain
 * even lower level structures, but in logical level it should be
 * atomic.
 *
 * The class should be abstract, but is not because @ref Array doesn't
 * allow that.
 **/
class Gene : public Genstruct {
	decl_dynamic (Gene);
  public:

	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
					Gene		() {FORBIDDEN}
					Gene		(const GeneticID& n, double mut=1);
					Gene		(const Gene& o) : Genstruct (o) {shallowCopy (o);}

	/** Sets the mutation coefficient */
	void						setMutability	(double val) {mutability=val;}

	/** Returns the mutation coefficient. */
	double						getMutability	() const {return mutability;}

	// Implementations
		
	/** Implementation for @ref Genstruct */
	virtual void				addPrivateGenes	(Gentainer& g, const StringMap& params) {;}

	/** Implementation for @ref Genstruct */
	virtual void				copy			(const Genstruct& other);

	/** Implementation for @ref Genstruct */
	virtual const Genstruct*	getGene			(const GeneticID& nam) const {
		return (id==nam)? this : (const Gene*) NULL;
	}
	virtual void				check			() const {Genstruct::check ();}
	
  protected:
	/** A coefficient for the mutation propability OR distribution
	 * function. Its semantics are dependent on the specific gene
	 * subclass.
	 **/
	double	mutability;

	virtual int					calc_len	() const {return 1;}
	void						shallowCopy	(const Gene& o) {mutability=o.mutability;}

  private:
	Gene& operator= (const Gene& orig) {FORBIDDEN; return *this;}
};

// The Gene is your friend!

/** Macro to make making certain gene subclasses easier (see the header).
 **/
#define DefGene(gclass) \
class gclass : public Gene { \
	decl_dynamic (gclass);\
  public:\
						gclass				(const GeneticID& name=NULL) : Gene (name) {;}\
						gclass				(const gclass& o) : Gene (o) {;}\
	void				init				() {;}\
	bool				pointMutate			(const MutationRate& k) {return false;}\
	void				copy				(const Genstruct& o) {Gene::copy(o);}\
	Genstruct*			replicate			() const {return new gclass (id);}\
	bool				execute				(const GeneticMsg& msg) const;\
	void				addPrivateGenes		(Gentainer& g, const StringMap& params);\
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           ----  o                        ----                             //
//           |   )     _    ___            |      ___    _    ___            //
//           |---  | |/ \   ___| |/\ \   | | --- /   ) |/ \  /   )           //
//           |   ) | |   | (   | |    \  | |   \ |---  |   | |---            //
//           |___  | |   |  \__| |     \_/ |___/  \__  |   |  \__            //
//                                    \_/                                    //
///////////////////////////////////////////////////////////////////////////////

/** A binary gene with values 0 and 1.
 **/
class BinaryGene : public Gene {
	decl_dynamic (BinaryGene);
  public:

	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
							BinaryGene	() {FORBIDDEN}

	/** The standard constructor.
	 *
	 *  @param id Name of the gene.
	 *  @param mut Local mutation rate coefficient.
	 *
	 *  @param initP Probability of value 1 in initialization. By
	 *  default the probability of 0 and 1 are uniform (0.5). You can,
	 *  of course, use this to set the initial value: use 0.0 for
	 *  value 0, and 1.0 for value 1.
	 **/
							BinaryGene	(const GeneticID& id,
										 double mut=1.0, double initP=0.5);

							BinaryGene	(const BinaryGene& o) : Gene (o) {shallowCopy (o);}

	// Access operators
	
	/** Sets the value of the gene. */
	BinaryGene&				set			(bool value) {mValue = value; return *this;}

	/** Returns the value of the gene as either true (1) or false (0). */
	bool					getvalue	() const {return mValue;}

	/** Sets the probability of value 1 in initialization. The
	 *  default initialization probability (see the constructor) is 0.5.
	 **/
	void					setInitP	(double initP) {mInitP=initP;}

	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual void			init		();

	/** Implementation for @ref Genstruct.
	 *
	 *  Flips the current value over with propability k*c, where c is
	 *  the gene-local probability coefficient, and k is described
	 *  below.
	 *
	 *  @param k Mutation rate parameters. Only the binaryrate
	 *  parameter is usefor for us.
	 **/
	virtual bool			pointMutate	(const MutationRate& k);
	/** Implementation for @ref Genstruct. The distance is calculated
	 *  as trivial case of Hamming distance.
	 **/
	virtual double			equality	(const Genstruct& other) const {
		return (mValue == static_cast<const BinaryGene&> (other).mValue)? 1.0 : 0.0;
	}
	/** Implementation for @ref Genstruct */
	virtual void			copy		(const Genstruct& other);
	/** Implementation for @ref Genstruct */
	virtual Genstruct*		replicate	() const {return new BinaryGene (*this);}
	/** Implementation for @ref Genstruct */
	virtual void			print		(TextOStream& out) const;
	/** Implementation for @ref Object */
	virtual DataOStream&	operator>>	(DataOStream& out) const;

  private:
	void					shallowCopy	(const BinaryGene& o) {
		mValue=o.mValue;
		mInitP=o.mInitP;
	}
	BinaryGene& operator= (const BinaryGene& orig) {FORBIDDEN; return *this;}

	bool	mValue;
	double	mInitP;
};

// A dead binary is the best binary



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//       _               ----- |                 ----                       //
//      / \    _         |     |       ___   |  |      ___    _    ___      //
//     /   \ |/ \  \   | |---  |  __   ___| -+- | --- /   ) |/ \  /   )     //
//     |---| |   |  \  | |     | /  \ (   |  |  |   \ |---  |   | |---      //
//     |   | |   |   \_/ |     | \__/  \__|   \ |___/  \__  |   |  \__      //
//                  \_/                                                     //
//////////////////////////////////////////////////////////////////////////////

/** Abstract interface class for genes that encode a floating-point
 *  value.
 *
 *  You can read about the creation of floating-point genes and how
 *  they are inserted in the genomes from the specific subclass
 *  documentation. The values of the genes are read from the genome as
 *  follows:
 *
 *  double x = dynamic_cast<const AnyFloatGene&> (*individual.getGene("myGene")).getvalue();
 *
 *  This class is actually not abstract, because of technical
 *  limitations in the dynamic runtime information system.
**/
class AnyFloatGene : public Gene {
	decl_dynamic (AnyFloatGene);
  public:
	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
					AnyFloatGene	() {FORBIDDEN}
	
	/** The subclass should call this constructor to set the name, value range and
	 * mutation parameter coefficient.
	 *
	 * @param id Name of the gene.
	 * @param min Lower limit for the floating-point value range.
	 * @param max Upper limit for the floating-point value range.
	 *
	 * @param m Mutation parameter coefficient; the exact semantics of
	 * this variable depends on the subclass.
	 **/
					AnyFloatGene	(const GeneticID& id, double min, double max, double m=1.0);
					AnyFloatGene	(const AnyFloatGene& o) : Gene (o) {
						shallowCopy (o);
					}

	/**
	 * Returns the decoded, phenotypic value of the gene
	 *
	 * The inheritor MUST overload this method.
	 *
	 * @exception must_overload
    **/
	virtual double	getvalue	() const {MUST_OVERLOAD; return 0.0;}
	
	// Implementations

	/** Standard copy operator. Implementation for @ref Genstruct. */
	virtual void	copy		(const Genstruct& other);

	/** Standard check operator. Implementation for @ref Object.
	 *  @exception assertion_failed check_failed
	 **/
	virtual void	check		() const;
	
  protected:
	/** The lower limit for the gene value */
	double			mMin;

	/** The upper limit for the gene value */
	double			mMax;
	
  private:
	/** Forbid the copy operator */
	AnyFloatGene&	operator= (const AnyFloatGene& orig) {FORBIDDEN; return *this;}
	
	void			shallowCopy	(const AnyFloatGene& o);

};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               ----- |                 ----                               //
//               |     |       ___   |  |      ___    _    ___              //
//               |---  |  __   ___| -+- | --- /   ) |/ \  /   )             //
//               |     | /  \ (   |  |  |   \ |---  |   | |---              //
//               |     | \__/  \__|   \ |___/  \__  |   |  \__              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// External
class FloatMutator;

/** Floating-point gene with native floating point representation;
 *  implementation of Evolution Strategies real-valued genes.
 *
 *  The mutation operation mutates the value every time by a small,
 *  normally distributed value (that can be global or self-regulated).
**/
class FloatGene : public AnyFloatGene {
  public:

	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
							FloatGene	() {FORBIDDEN}

	/** Standard constructor. Usually used for example:
	 *
	 * myGenome.add (new FloatGene ("myGene", 0.2, 0.8, 0.1));
	 *
	 * @param id Name of the gene.
	 * @param min Lower limit for the floating-point value range.
	 * @param max Upper limit for the floating-point value range.
	 * @param m Mutation coefficient.
	 **/
							FloatGene	(const GeneticID& id, double min, double max,
										 double m=1.0, double value=-999);
							FloatGene	(const FloatGene& o) : AnyFloatGene (o) {
								shallowCopy (o);
							}

	/** Sets the gene value. */
	FloatGene&				set			(double val) {value = val; return *this;}
	
	/** Implementation for @ref AnyFloatGene. Returns the "decoded"
	 *  value of the gene. There actually isn't any "decoding" in this
	 *  gene class, and the phenotypic value equals the genotypic.
	 **/
	virtual double			getvalue	() const {return value;}
	

	/** Sets the value range to be handled as circular; if a value
	 *  near the upper limit is mutated, the result can be a small
	 *  value.
	 *
	 * @return Self.
	 **/
	FloatGene&				setCircular	(bool c) {mCircular=c; return *this;}

	/** Sets a customized mutation operator for the gene. See
	 *  inheritors of @ref FloatMutator for additional
	 *  information. The default for @ref FloatGene is to use gaussian mutation.
	 *
	 * @return Self.
	 **/
	FloatGene&				setMutator	(FloatMutator* m) {mMutator=m; return *this;}

	// Implementations

	/** Implementation for @ref Genstruct.  The operator initializes
	 *  the gene with a random variable with uniform distribution.
	 **/
	virtual void			init		();

	/** Implementation for @ref Genstruct.
	 *
	 *  We use two control parameters. First is the mutation rate,
	 *  which is a probability to make the actual mutation. This
	 *  parameter is our own addition; with Evolution Strategies the
	 *  mutation rate is always 1.0 so that a mutation always occurs.
	 *
	 *  The other control parameter is mutation variance, which
	 *  controls the usually normally distributed mutation. If,
	 *  however, the mutation operator has been redefined with
	 *  setMutator(), the variance-parameter may be used with
	 *  different semantics.
	 *
	 *  @param k Mutation rate coefficients; mutation rate and variance.
	 **/
	virtual bool			pointMutate	(const MutationRate& k);

	/** Implementation for @ref Genstruct.
	 *
	 *  Genetic distance is measured as real-valued distance between
	 *  the genes. Note that, if the distance is measured for a set of @ref
	 *  FloatGene genes, the distance will be orthogonal, not
	 *  euclidean.
	 **/
	virtual double			equality	(const Genstruct& other) const;
	virtual void			copy		(const Genstruct& other);
	virtual Genstruct*		replicate	() const {return new FloatGene (*this);}
	virtual void			print		(TextOStream& out) const;
	virtual void			check		() const;
	
  protected:

	/** The genotypic value of the gene */
	double			value;

	/** Local variance modifier for mutation variance. */
	double			mVariance;

	/** Should the mutations be done in circular domain? */
	bool			mCircular;

	/** Mutation method */
	FloatMutator*	mMutator;

  private:
	void					shallowCopy	(const FloatGene& o);
	FloatGene&				operator=	(const FloatGene& orig) {FORBIDDEN; return *this;}
	decl_dynamic (FloatGene);
};

// Gee, genes don't sink; they really do double



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ----  o     ----- |                 ----                          //
//        |   )    |  |     |       ___   |  |      ___    _    ___         //
//        |---  | -+- |---  |  __   ___| -+- | --- /   ) |/ \  /   )        //
//        |   ) |  |  |     | /  \ (   |  |  |   \ |---  |   | |---         //
//        |___  |   \ |     | \__/  \__|   \ |___/  \__  |   |  \__         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** Floating point value encoded genetically as a binary vector. This
 *  is the way to encode "real" values in canonical Genetic
 *  Algorithms.
 *
 *  The encoding has two modes: linear and Gray-encoded. In linear
 *  mode the value is simply a sum of bits: SUM(i=0..n-1) b(i)*2^i,
 *  where b(i) is the i:th bit in the genetic sequence b.
 **/
class BitFloatGene : public AnyFloatGene {
  public:
	
	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
							BitFloatGene	() {FORBIDDEN}

	/** Standard constructor. Usually used for example:
	 *
	 *  myGenome.add (new BitFloatGene ("myGene", 0.2, 0.8, 10, params, 1.0));
	 *
	 *  @param id Name of the gene.
	 *  @param min Lower limit for the floating-point value range.
	 *  @param max Upper limit for the floating-point value range.
	 *  @param bits Number of bits for encoding the value.
	 *  @param params Additional dynamic parameters.
	 *  @param params["grayCoding"] Should Gray coding be used. [Default:0]
	 *  @param m Mutation coefficient. [Default:1.0]
	 **/
							BitFloatGene	(const GeneticID& id, double min, double max,
											 int bits, const StringMap& params, double m=1.0);
							BitFloatGene	(const BitFloatGene& o) : AnyFloatGene (o) {
								shallowCopy (o);
							}

	virtual double			getvalue		() const;
	BitFloatGene*			setGrayCoding	(bool isg=true) {mGrayCoded=isg; return NULL;}
	
	// Implementations
		
	/** Initializes the bits totally randomly. */
	virtual void			init		() {mBits.init();}
	
	/** The point mutation forwards the mutation to the @ref Gentainer
	 *  holding the @ref BinaryGene bits that encode the integer value.
	 **/
	virtual bool			pointMutate	(const MutationRate& k) {return mBits.pointMutate (k);}

	/** The genetic distance is calculated as the Hamming distance
	 *  between the bits that encode the integer value. This might
	 *  really not be what you usually want, as some bits are usually
	 *  more important than others. It might be more appropriate to
	 *  use the difference between getvalue() decoded values.
	 **/
	virtual double			equality	(const Genstruct& o) const {mBits.equality(o); return 0.0;}
	virtual void			copy		(const Genstruct& other);
	virtual Genstruct*		replicate	() const {return new BitFloatGene (*this);}
	virtual void			print		(TextOStream& out) const;
	virtual void			recombine	(const Genstruct& a, const Genstruct& b);
	virtual void			check		() const;

  protected:
	/** Actual @ref BinaryGene sequence encoding the floating-point value.
	 **/
	Gentainer		mBits;

	/** Number of bits in the genetic representation of the
	 *  floating-point value. This could be calculated from the mBits
	 *  @ref Gentainer, but that would be slow, so we cache the value
	 *  here.
	 **/
	int				mBitCount;

	/** Gray coding mode flag. */
	bool			mGrayCoded;
	
	virtual int				calc_len	() const {return ((const Genstruct&)mBits).calc_len();}
	
  private:
	void					shallowCopy	(const BitFloatGene& o);
	BitFloatGene&			operator=	(const BitFloatGene& orig) {FORBIDDEN; return *this;}
	decl_dynamic (BitFloatGene);
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//            _               ---            ----                            //
//           / \    _          |    _    |  |      ___    _    ___           //
//          /   \ |/ \  \   |  |  |/ \  -+- | --- /   ) |/ \  /   )          //
//          |---| |   |  \  |  |  |   |  |  |   \ |---  |   | |---           //
//          |   | |   |   \_/ _|_ |   |   \ |___/  \__  |   |  \__           //
//                       \_/                                                 //
///////////////////////////////////////////////////////////////////////////////

/** Abstract interface for any genes that encode an integer value.
 *
 *  You can read about the creation of integer-valued genes and how
 *  they are inserted in the genomes from the specific subclass
 *  documentation. The values of the genes are read from the genome as
 *  follows:
 *
 *  double x = dynamic_cast<const AnyIntGene&> (*individual.getGene("myGene")).getvalue();
 *
 *  Because of technical limitations, the class is not actually
 *  abstract.
 **/
class AnyIntGene : public Gene {
	decl_dynamic (AnyIntGene);
  public:

	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
					AnyIntGene	() {FORBIDDEN}
	
	 /** The subclass should call this constructor to set the name and
	 * some variables.
	 *
	 * @param id Name of the gene.
	 * @param min Lower limit for the floating-point value range.
	 * @param max Upper limit for the floating-point value range.
	 * @param m Mutation coefficient.
	 **/
					AnyIntGene	(const GeneticID& id, int min, int max, double m=1.0);
	
					AnyIntGene	(const AnyIntGene& o) : Gene (o) {
						shallowCopy (o);
					}

	/**
	 * Returns the decoded phenotypic value of the gene.
	 *
	 * The inheritor MUST overload this method.
    **/
	virtual int		getvalue	() const {MUST_OVERLOAD; return 0;}
	
	// Implementations

	/** Standard copy operator. Implementation for @ref Genstruct. */
	virtual void	copy		(const Genstruct& other);
	
	/** Standard check operator. Implementation for @ref Object.
	 *  @exception assertion_failed check_failed
	 **/
	virtual void	check		() const {ASSERT(mMin<mMax); Gene::check ();}
	
  protected:
	/** The lower limit for the gene value */
	int				mMin;

	/** The upper limit for the gene value */
	int				mMax;

	void			shallowCopy	(const AnyIntGene& o);
};


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                    ---            ----                                    //
//                     |    _    |  |      ___    _    ___                   //
//                     |  |/ \  -+- | --- /   ) |/ \  /   )                  //
//                     |  |   |  |  |   \ |---  |   | |---                   //
//                    _|_ |   |   \ |___/  \__  |   |  \__                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** Integer (or multivalue) gene with independent values.
 *
 * This is the integer equivalent of the @ref FloatGene. The gaussian mutation
 * method may not suit this data type so well, but... that's only life.
**/
class IntGene : public AnyIntGene {
	decl_dynamic (IntGene);
	int		mValue;
  public:

	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
							IntGene		() {FORBIDDEN}
	
							IntGene		(const GeneticID& id, int min, int max,
										 double mut=1,int value=-9999);
							IntGene		(const IntGene& o) : AnyIntGene (o) {shallowCopy (o);}

	/** Sets the value of the gene. */
	virtual IntGene&		set			(int val) {mValue = val; return *this;}

	/** Implementation for @ref AnyIntGene.
	 *
	 * @return The decoded phenotypic value of the gene.
	 **/
	virtual int				getvalue	() const;

	// Implementations

	/** This implementation initializes the gene into any value with
	 *  equal probability.
	 **/
	virtual void			init		();

	/** Implementation for @ref Genstruct. If the mutation actualizes,
	 *  it mutates the integer gene into any other value with uniform
	 *  probability.
	 *
	 * @return If a mutation really actualized or not.
	 **/
	virtual bool			pointMutate	(const MutationRate& k);

	/** Implementation for @ref Genstruct. The distance measure for
	 *  two integer genes is 0.0 if they are equal, and 1.0 if they
	 *  are different.
	 **/
	virtual double			equality	(const Genstruct& other) const;
	virtual void			copy		(const Genstruct& other);
	virtual Genstruct*		replicate	() const {return new IntGene (*this);}
	virtual void			print		(TextOStream& out) const;
	virtual void			check		() const;

  private:
	void					shallowCopy	(const IntGene& o);
	IntGene&				operator=	(const IntGene& orig) {FORBIDDEN; return *this;}
};

// IntGene? .. hmm .. Wouldn't InGen be better?


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//             ----  o     ---            ----                               //
//             |   )    |   |    _    |  |      ___    _    ___              //
//             |---  | -+-  |  |/ \  -+- | --- /   ) |/ \  /   )             //
//             |   ) |  |   |  |   |  |  |   \ |---  |   | |---              //
//             |___  |   \ _|_ |   |   \ |___/  \__  |   |  \__              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** Integer gene that is implemented as a vector of @ref BinaryGene
 *  genes.
 *
 *  Because of the vector representation, the values are correlated,
 *  which would ne be case with @ref IntGene.
 *
 *  The range of the gene values must be a power of two. For example
 *  range [0..1] is 2=2^1, and range [10..265] is 256=2^8.
 **/
class BitIntGene : public AnyIntGene {
  public:
	/** Default constructor, FORBIDDEN! Exists only because of RTTI system.
	 **/
							BitIntGene	() {FORBIDDEN}

	/** Standard way to construct the gene.
	 *
	 * @param id Name of the gene.
	 * @param min Minimum value.
	 * @param max Maximum value.
	 * @param bits Number of bits in the genetic representation. This must equal sqrt(man-min).
	 * @param params Additional dynamic parameters.
	 * @param params["grayCoding"] Should Gray coding be used for the value? [Default:0]
	 * @param m Mutation probability coefficient. [Default:1.0]
	 **/
							BitIntGene	(const GeneticID& id, int min, int max,
											 int bits, const StringMap& params, double m=1.0);
							BitIntGene	(const BitIntGene& o) : AnyIntGene (o) {
								shallowCopy (o);
							}

	/** Returns the decoded phenotypic value of the gene.
	 **/
	virtual int				getvalue		() const;

	/** Sets the Gray coding mode on (true) or off (false).
	 *
	 * The return value makes it possible to use the method as
	 * follows:
	 *
	 * myGenome.add ((new BitIntGene ("myGene", 0, 7, 3, params))->setGrayCoding (true));
	 *
	 * @return Returns self.
	 **/
	BitIntGene*				setGrayCoding	(bool gray=true) {mGrayCoded=gray; return this;}
	
	// Implementations

	/** Initializes the bits totally randomly. */
	virtual void			init		() {mBits.init();}

	/** The point mutation forwards the mutation to the @ref Gentainer
	 *  holding the @ref BinaryGene bits that encode the integer value.
	 **/
	virtual bool			pointMutate	(const MutationRate& k) {return mBits.pointMutate (k);}

	/** The genetic distance is calculated as the Hamming distance
	 *  between the bits that encode the integer value. If the values
	 *  of the bits do not have an equal significanse for the
	 *  phenotype, this way of calculating may not be what you want.
	 **/
	virtual double			equality	(const Genstruct& o) const {return mBits.equality(o);}
	virtual void			copy		(const Genstruct& other);
	virtual Genstruct*		replicate	() const {return new BitIntGene (*this);}
	virtual void			print		(TextOStream& out) const;
	virtual void			recombine	(const Genstruct& a, const Genstruct& b);
	virtual void			check		() const;

  protected:

	/** Binary gene vector. */
	Gentainer	mBits;

	/** Number of bits in the vector. This can be calculated from the
	 *  mBits @ref Gentainer, but that would be slow, so we cache the
	 *  value here.
	 **/
	int			mBitCount;

	/** Gray coding mode flag. */
	bool		mGrayCoded;
	
	virtual int				calc_len	() const {return ((const Genstruct&)mBits).calc_len();}
	
  private:
	void					shallowCopy	(const BitIntGene& o);
	BitIntGene&				operator= (const BitIntGene& orig) {FORBIDDEN; return *this;}
 	decl_dynamic (BitIntGene);
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//              ---                      ----                                //
//               |    _    |   ___      |      ___    _    ___               //
//               |  |/ \  -+- /   ) |/\ | --- /   ) |/ \  /   )              //
//               |  |   |  |  |---  |   |   \ |---  |   | |---               //
//              _|_ |   |   \  \__  |   |___/  \__  |   |  \__               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** A reference gene. Executing this gene automatically executes
 *  another gene.
 **/
class InterGene : public Gene {
	decl_dynamic (InterGene);
	GeneticID	targetGene;
  public:
	/** Standard constructor.
	 *
	 *  @param name Name of the gene.
	 *  @param name trg Name of the target gene which this gene acts as a reference to.
	 **/
							InterGene	(const GeneticID& name=NULL, const GeneticID& trg=NULL)
									: Gene (name), targetGene (trg) {;}
							InterGene	(const InterGene& o) : Gene (o) {shallowCopy (o);}

	// Implementations

	/** Implementation for @ref Genstruct. The initialization command
	 *  is not passed through to the referred gene.
	 **/
	virtual void			init		() {;}
	virtual Genstruct*		replicate	() const {return new InterGene (*this);}
	virtual void			copy		(const Genstruct& other);
	virtual bool			execute		(const GeneticMsg& msg) const;
	virtual void			print		(TextOStream& out) const;
	virtual bool			pointMutate	(const MutationRate& k) {return false;}

  private:
	void					shallowCopy	(const InterGene& o) {targetGene=o.targetGene;}
	InterGene&				operator=	(const InterGene& orig) {FORBIDDEN; return *this;}
};



// I don't want any f*cking genes in my food.

#endif
