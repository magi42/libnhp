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

/*
  Neljänneksi tulee Tuulen kirja. Tämä kirja ei ole yhteydessä minun
  Ichi-kouluuni vaan muihin strategian kouluihin. Tuulella tarkoitan
  strategian vanhoja perinteitä, tämän päivän perinteitä ja
  perhetraditioita. Siten voi selkeästi selittää maailman
  strategiat. Tämä on perinne. On vaikeaa tuntea itseäsi jos et tunne
  muita. Kaikilla Teillä on sivupolkuja. Jos opiskelet Tietä
  päivittäin ja henkesi joutuu harhaan, voit luulla noudattavasi hyvää
  Tietä vaikka objektiivisesti se ei ole oikea Tie. Jos seuraat oikeaa
  Tietä mutta harhaannut hieman, tulee siitä myöhemmin suuri
  poikkeama. Sinun täytyy tajuta tämä. Muut strategiat ovat tulleet
  otettaviksi pelkkänä miekkailuna, eikä se ole väärin että näin
  olisi. Minun strategiani hyöty, vaikkakin se sisältää miekkailun,
  sijaitsee erillisessä periaatteessa. Olen selittänyt mitä
  tavallisesti tarkoitetaan strategialla muissa kouluissa Perinteen
  (Tuulen) kirjassa.  */

#ifndef __GENETICS_H__
#define __GENETICS_H__

#include <magic/mobject.h>
#include <magic/mstring.h>
#include <magic/mmap.h>

using namespace MagiC;

// Externals
// class MagiC::DataOStream;
// class MagiC::OStream;
class Individual;
class Gene;
class InterGene;

// Internals
class Genstruct;
class Gentainer;

enum printflags {PRINT_CLASSNAMES=0, PRINT_NAMES, PRINT_HIDDEN};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//            ----                       o       |   |                       //
//           |      ___    _    ___   |     ___  |\ /|  ____                 //
//           | --- /   ) |/ \  /   ) -+- | |   \ | V | (      ___            //
//           |   \ |---  |   | |---   |  | |     | | |  \__  (   \           //
//           |___/  \__  |   |  \__    \ |  \__/ |   | ____)  ---/           //
//                                                            __/            //
///////////////////////////////////////////////////////////////////////////////

typedef String GeneticID;

/** A genetic message that activates genes. 
 **/
class GeneticMsg {
  public:

	/** Create a genetic message.
	 *
	 * @param rcvr Receiver of the message. Typically a gene name as a
	 * @ref String. It is planned that the @ref GeneticID could also
	 * be a name mask, and also perhaps something else than a @ref
	 * String.
	 *
	 * @param ind Individual. See member variable "host" below.
	 **/
	GeneticMsg	(const GeneticID& rcvr, Individual& ind) : mrReceiver (rcvr), mrHost (ind) {
		mpPrevious = NULL;
	}

	/** Copy constructor. */
	GeneticMsg (const GeneticMsg& orig) :
			mrReceiver (orig.mrReceiver), mrHost (orig.mrHost) {}

	virtual ~GeneticMsg () {
	}

	/** Pushes the current caller into call stack.
	 *
	 * Use this like as follows:
	 *
	 * gmsg.pushCaller (this);
	 *
	 * child.execute (gmsg);
	 *
	 * Since the caller is not currently implemented with stack, it should
	 * be repushed if another children are successively executed.
	 *
	 * (This method is virtual only because we want to have a virtual table for the class!)
	 **/
	virtual void	pushCaller	(Genstruct* pClr) {mpPrevious = pClr;}

	/** Returns the caller that was previously pushed */
	Genstruct*	getCaller	() const {return mpPrevious;}

	/** The name of the gene that is supposed to receive this message. */
	const GeneticID&	mrReceiver;

	/** The individual whose genome we are handling.
	 *
	 * This information is available for genes to allow them to decode
	 * themselves as phenotypes of the individual. It is also possible
	 * to activate or seek other genes of the host.
	 **/
	Individual&			mrHost;

  private:
	/** Refers to the previous handler of this message. We are lazy
	 * and don't want to use stack.
	 **/
	Genstruct*	mpPrevious;

	GeneticMsg& operator= (const GeneticMsg& orig) {FORBIDDEN; return *this;}

	friend class InterGene;
	friend class Individual;
};

/** Hmm, could this be active? Hmm, could this be genetically encoded?
	Hmm, many deep questions. **/


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//       |   |                     o            ----                        //
//       |\ /|        |   ___   |           _   |   )  ___   |   ___        //
//       | V | |   | -+-  ___| -+- |  __  |/ \  |---   ___| -+- /   )       //
//       | | | |   |  |  (   |  |  | /  \ |   | | \   (   |  |  |---        //
//       |   |  \__!   \  \__|   \ | \__/ |   | |  \   \__|   \  \__        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** A parameter-passing class for passing mutation rates to
 * atomic genes @ref BinaryGene, @ref AnyIntGene, and @ref AnyFloatGene.
 **/
class MutationRate {
	double			mBinaryRate;
	double			mIntRate;
	double			mFloatRate;
	double			mFloatVariance;
	bool			mOneBitMutation;
	bool			mAutoAdaptation;
	//FloatMutator*	mFloatMutator;
  public:
						MutationRate	() {
							mBinaryRate=mIntRate=mFloatRate=mFloatVariance=0.0;
							mOneBitMutation=false;
							mAutoAdaptation=false;
						}

						MutationRate	(const MutationRate& o) {
							// Do copying quickly
							memcpy (this, &o, sizeof (MutationRate));
						}

	/** Constructs the object by multiplying two existing mutation
	 * rates.
	 **/
						MutationRate	(const MutationRate& o, const MutationRate& m) {
							mBinaryRate = o.mBinaryRate * m.mBinaryRate;
							mIntRate = o.mIntRate * m.mIntRate;
							mFloatRate = o.mFloatRate * m.mFloatRate;
							mFloatVariance = o.mFloatVariance * m.mFloatVariance;
							mOneBitMutation = o.mOneBitMutation || m.mOneBitMutation;
							mAutoAdaptation = o.mAutoAdaptation || m.mAutoAdaptation;
						}
	
	/** Reads the mutation rates from the given gentainer
	 **/
						MutationRate	(const Gentainer& g);
	
	/** Adds the genes controlling the mutation rates to the given gentainer. */
	static void			addGenesTo		(Gentainer& g, const StringMap& params);

	/** Returns the mutation rate for @ref BinaryGene. */
	double				binaryRate		() const {return mBinaryRate;}
	/** Returns the mutation rate for @ref IntGene. */
	double				intRate			() const {return mIntRate;}
	/** Returns the mutation rate for @ref FloatGene. */
	double				doubleRate		() const {return mFloatRate;}
	/** Returns the mutation variance for @ref FloatGene. */
	double				doubleVariance	() const {return mFloatVariance;}

	/** Sets the mutation rate for @ref BinaryGene */
	void				binaryRate		(double x) {mBinaryRate=x;}
	/** Sets the mutation rate for @ref IntGene */
	void				intRate			(double x) {mIntRate=x;}
	/** Sets the mutation rate for @ref FloatGene */
	void				doubleRate		(double x) {mFloatRate=x;}
	/** Sets the mutation variance for @ref FloatGene */
	void				doubleVariance	(double x) {mFloatVariance=x;}

	/** Order that only a SINGLE leaf-level genstruct should be mutated.
	 **/
	void				oneBitMutation	(bool mt) {mOneBitMutation = mt;}

	/**
	 * @return The value of the oneBitMutation flag.
	 **/
	bool				oneBitMutation	() const {return mOneBitMutation;}

	/** Setting this to TRUE causes the message-passing to read the
        autoadaptation genes from the @ref Gentainer objects.
	**/
	void				autoAdaptation	(bool aa) {mAutoAdaptation=aa;}

	/** Returns the state of autoadaptivity */
	bool				autoAdaptation	() const {return mAutoAdaptation;}

  private:
	MutationRate& operator= (const MutationRate& orig) {FORBIDDEN; return *this;}
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//               ----                                                        //
//              |      ___    _    ____  |             ___   |               //
//              | --- /   ) |/ \  (     -+- |/\ |   | |   \ -+-              //
//              |   \ |---  |   |  \__   |  |   |   | |      |               //
//              |___/  \__  |   | ____)   \ |    \__!  \__/   \              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

EXCEPTIONCLASS (invalid_gene_name);

/** The highest genetic abstraction: a genetic structure. These
* structures can be "genes" or "chromosomes" or "genomes", or whatever.
**/
class Genstruct : public Object {
	decl_dynamic (Genstruct);
  public:

	/** Creates the structure with the given name. The name can be
	 * used to access the structure from @ref Gentainer containers.
	 **/
								Genstruct	(const GeneticID& name=NULL);
								Genstruct	(const Genstruct& other);


	/** Gene access operator.
	 * @return The gene with the given name within the structure.
	 * @throw invalid_gene_name if the gene was not found within the structure.
	 **/
	const Genstruct&			operator[]	(const GeneticID& name) const {
		return *getGene (name);
	}

	/** Gene access operator, just like above, but for non-const objects.
	 * @return The gene with the given name within the structure.
	 **/
	Genstruct&					operator[]	(const GeneticID& name) {
		return const_cast<Genstruct&>(const_cast<const Genstruct*>(this)->operator[] (name));
	}

	/** Retrieves the name (ID) of the structure */
	const GeneticID&			getID		() const {return id;}

	/** Calculates recursively the total number of genes within the
     *  genetic structure.
	**/
	int							length		() const;

	//////////////////////////////////////////////////////////////////////
	// Virtual functions
	
	/** Generates an initial, random state for the structure.
	 **/
	virtual void				init		() {MUST_OVERLOAD}

	/** Fetches a gene from the structure by it's name.
	 *
	 * To be logical, this really should be a method of @ref
	 * Gentainer, but there are cases when some classes can't inherit
	 * @ref Gentainer, and they still need this.
	**/
	virtual const Genstruct*	getGene		(const GeneticID& name) const {MUST_OVERLOAD; return NULL;}

	/** Transmits a sad order to execute a gene. The receiver of the
	 * message reads in the upper left corner of the message. It
	 * should then execute itself. It's good that objects don't have
	 * feelings. (Hmm, they don't?)
	 **/
	virtual bool				execute		(const GeneticMsg& msg) const {return false;}

	/** Adds internal genes to itself.
	 *
	 * @param g Owner @ref Gentainer.
	 **/
	virtual void				addPrivateGenes	(Gentainer& g) {MUST_OVERLOAD}
	/** Adds internal genes to itself. Additional parameters are
	 * passed in the @ref String @ref Map.
	 *
	 * @param g Owner @ref Gentainer.
	 * @param params Dynamic parameters in a @ref String @ref Map.
	 **/
	virtual void				addPrivateGenes	(Gentainer& g, const StringMap& params) {MUST_OVERLOAD}
	
	/** Mutates the structure.
	 *
	 * @param r Mutation rate coefficients.
	 *
	 * @return TRUE if an actual mutation (change on genotype) has
	 * occurred.
	 **/
	virtual bool				pointMutate	(const MutationRate& r) {MUST_OVERLOAD; return false;}
	
	/** Makes this structure a recombination of given parent
	 * structures. If no internal recombination actualizes within the
	 * structure, the parent 'a' must always be copied.
	 **/
	virtual void				recombine	(const Genstruct& a, const Genstruct& b) {copy(a);}

	/** Compares how genotypically similar the structure is to
	 * another. The units of the distance are usually
	 * application-dependent. One might, for example, use same
	 * euclidean distance or Hamming distance measurement with all
	 * genes, but then some genes may really be irrelevant for some
	 * problem, and they may have to be weighted. Thus, implementing
	 * this method in a generic way may be rather difficult.
	 *
	 * @return The distance between the structures.
	 **/
	virtual double				equality	(const Genstruct& other) const {MUST_OVERLOAD; return 0.0;}

	/** Clones a copy of self, RECURSIVELY. It's funny how simple reproduction
	 * is nowadays.
	 **/
	virtual Genstruct*			replicate	() const {MUST_OVERLOAD; return NULL;}

	/** Copies another structure to self (recursively). The other
	 * structure must be structurally EXACTLY equivalent to self
	 * (for example, replicated previously). This is cool.
	 **/
	virtual void				copy		(const Genstruct& other) {MUST_OVERLOAD}

	/** Recursively prints the genome to the given stream. This is
	 * most cool.
	 **/
	virtual void				print		(TextOStream& out) const;

	// The world has become so virtual that it's a wonder we any more
	// know what is real. For example, this class is abstract.

	/** Another structural output operator, implementation from @ref
	 * Object.
	 **/
	virtual DataOStream&		operator>>	(DataOStream& out) const;

	/** Sets the hide-flag of the structure, for example for printing
	 * methods, etc.
	 **/
	Genstruct&					hide		(bool h=true) {mHidden=h; return *this;}

	/** Is the genstruct hidden from printing? */
	bool						isHidden	() const {return mHidden;}

	/** Actually calculates the true length of the genome (recursively).
	 **/
	virtual int					calc_len	() const {MUST_OVERLOAD; return 0;}

	// Implementations

	/** Implementation for @ref Object. */
	virtual void				check		() const;
	
  protected:

	/** Inheritors should use this method to copy the contents.
	 **/
	void						copyGenstr	(const Genstruct& other) {
		id = other.id;
		size = other.size;
		mHidden = other.mHidden;
	}

	/** The identification label of the structure. It might be a string
	 * or it might be a number. But I am a genstruct, I am not a number.
	 **/
	GeneticID					id;

  private:	
	/** Cached number of genetic atoms (genes, codons, bases, bits or
	 * whatever). This is survival of the fattest.
	 **/
	mutable int					size;

	/** A flag: should the structure be hidden in dumps (default=false);
	 **/
	bool						mHidden;

	Genstruct operator= (const Genstruct& orig) {FORBIDDEN; return *this;}

	/** The Gentainer is your friend. Love the Gentainer. */
	friend class Gentainer;
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                ----                       o                               //
//               |      ___    _    |   ___      _    ___                    //
//               | --- /   ) |/ \  -+-  ___| | |/ \  /   ) |/\               //
//               |   \ |---  |   |  |  (   | | |   | |---  |                 //
//               |___/  \__  |   |   \  \__| | |   |  \__  |                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** A genetic container. It is a structurally middle-level @ref
 * Genstruct that contains atomic genes as well as other
 * containers. Many genetic operations are implemented here.
 **/
class Gentainer : public Genstruct {
	decl_dynamic (Gentainer);
  public:

						Gentainer	(const GeneticID& name = NULL);

						Gentainer	(const Gentainer& orig);

	/** Appends a new substructure to the structure. */
	void				add			(Genstruct* newstruct);
	
	/** Fetches the gene with the given ID within the
	 * structure. This method is not recursive.
	 **/
	const Genstruct&	operator[]	(const GeneticID& name) const {
		return *getGene (name);
	}
	/** Fetches the gene with the given ID within the
	 * structure. This method is not recursive.
	 *
	 * Non-const version of the above method.
	 **/
	Genstruct&			operator[]	(const GeneticID& name) {
		return const_cast<Genstruct&>(const_cast<const Gentainer*>(this)->operator[](name));
	}

	/** Fetches a substructure by it's index number. */
	const Genstruct&	operator[]	(int i) const {return substructs[i];}
	/** Fetches a substructure by it's index number. Non-const version.*/
	Genstruct&			operator[]	(int i) {return substructs[i];}

	/** Returns the number of immediate substructures contained. */
	int					size		() const {return substructs.size();}

	/** Enable self-adjusting of mutation rates within the
	 * container. When it is enabled, mutation-rate genes are inserted
	 * in the genome and used as coefficient for all mutation orders.
	 **/
	void				selfadjust	(bool val=true) {self_adjust=val;}

	/** Sets the local recombination rate coefficient. */
	void				recombRate	(double rate) {mRecombRate=rate;}

	// Implementations

	virtual void				init		();
	virtual void				addPrivateGenes	(Gentainer& g, const StringMap& params);
	virtual const Genstruct*	getGene		(const GeneticID& name) const;
	virtual bool				pointMutate	(const MutationRate& k);
	virtual void				recombine	(const Genstruct& a, const Genstruct& b);
	virtual double				equality	(const Genstruct& other) const;
	virtual Genstruct*			replicate	() const;
	virtual void				copy		(const Genstruct& other);
	virtual void				print		(TextOStream& out) const;
	virtual bool				execute		(const GeneticMsg& msg) const;

	virtual DataOStream&		operator>>	(DataOStream& out) const;
	virtual void				check		() const;
	
  protected:
	virtual int			calc_len	() const;

	/** Substructures. */
	Array<Genstruct>	substructs;

	/** Self-adjusting (autoadaptive) mutation rates. Default: false. */
	bool				self_adjust;

	/** Recombination rate coefficient. Default: 1 */
	double				mRecombRate;

  private:
	Gentainer& operator= (const Gentainer& orig) {FORBIDDEN; return *this;}
};

//
// I CALL AND I CALL!
// But in vain -
// no genstruct hears,
// no genstruct cares.
//
// So vast is the genome,
// so full of abstract gentainers
// that my messages are lost,
// in unoverridden virtual methods.
//


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      ----                                                //
//                     |      ___    _               ___                    //
//                     | --- /   ) |/ \   __  |/|/| /   )                   //
//                     |   \ |---  |   | /  \ | | | |---                    //
//                     |___/  \__  |   | \__/ | | |  \__                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** The highest level container of genetic information. Typically the
 * top level genstructs within the genome are recombined uniformly, to
 * represent the independet assorment of chromosomes in biological
 * genomes.
**/
class Genome : public Gentainer {
	decl_dynamic (Genome);
  protected:
	/** How many times a phenotype of this genome has been the King of
	 * the Kang. Or King of the Jungle.
	 **/
	int	kings;

  public:	

								Genome		();
								Genome		(const Genome& other);
	
	/** Adds one score for an phenotype of this genome having scored
	 * the kinghood. Useful only in elitist models.
	 **/
	void						addking		() {kings++;}

	/** Returns the number of kinghoods attained so far. Useful only
	 * in elitist models.
	 **/
	int							getkings	() const {return kings;}

	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual void				init		();
	/** Implementation for @ref Genstruct. */
	//virtual void				recombine	(const Genstruct& a, const Genstruct& b);
	/** Implementation for @ref Genstruct. */
	virtual void				print		(TextOStream& out) const;
	/** Implementation for @ref Genstruct. */
	virtual void				addPrivateGenes (Gentainer& g, const StringMap& pars);

	/** Implementation for @ref Object. */
	virtual DataOStream&		operator>>	(DataOStream& out) const;
	/** Implementation for @ref Object. */
	virtual void				check		() const;

  private:
	Genstruct& operator= (const Genstruct& orig) {FORBIDDEN; return *this;}
};




//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// |   |                      |      o     |  ----                              //
// |\ /|        _         --  |            | |      ___    _               ___  //
// | V |  __  |/ \   __  |  ) |  __  |  ---| | --- /   ) |/ \   __  |/|/| /   ) //
// | | | /  \ |   | /  \ |--  | /  \ | (   | |   \ |---  |   | /  \ | | | |---  //
// |   | \__/ |   | \__/ |    | \__/ |  ---| |___/  \__  |   | \__/ | | |  \__  //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

/** Genome with single set of chromosomes -- UNDER DEVELOPMENT.
 **/
class MonoploidGenome : public Genome {
  public:

				MonoploidGenome		();
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ___   o      |      o     |  ----                                   //
//      |  \     --  |            | |      ___    _               ___       //
//      |   | | |  ) |  __  |  ---| | --- /   ) |/ \   __  |/|/| /   )      //
//      |   | | |--  | /  \ | (   | |   \ |---  |   | /  \ | | | |---       //
//      |__/  | |    | \__/ |  ---| |___/  \__  |   | \__/ | | |  \__       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/** Genome with multiple sets of chromosomes -- UNDER DEVELOPMENT.
 **/
class DiploidGenome : public Genome {
	MonoploidGenome	genome_a, genome_b;
	
  public:

								DiploidGenome	();

	// Implementations

	virtual void		init		();
	virtual void		addPrivateGenes	(Gentainer& g, const StringMap& params);
	virtual int			length		() const;
	virtual const Gene&	operator[]	(const char* name) const  throw (invalid_gene_name);
	virtual bool		pointMutate	(const MutationRate& k);
	virtual void		recombine	(const Genstruct& a, const Genstruct& b);
	virtual double		equality	(const Genstruct& other) const;
	virtual Genstruct*	replicate	() const;
	virtual void		copy		(const Genstruct& other);
	virtual void		print		(TextOStream& out) const;
  protected:
	virtual int			calc_len	() const;

  private:
	DiploidGenome& operator= (const DiploidGenome& orig) {FORBIDDEN; return *this;}
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// -----                        ----                       o                 //
//   |        ___    _    ____ |      ___    _    |   ___      _    ___      //
//   |   |/\  ___| |/ \  (     | --- /   ) |/ \  -+-  ___| | |/ \  /   ) |/\ //
//   |   |   (   | |   |  \__  |   \ |---  |   |  |  (   | | |   | |---  |   //
//   |   |    \__| |   | ____) |___/  \__  |   |   \  \__| | |   |  \__  |   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/** An abstract transparent gentainer class. Not currently used.
 **/
/*
class TransGentainer {
  public:
  private:
	Gentainer*	gentainer;
};
*/

// And Void became to exist, together with the Light,
// The Mother Ocean was born, together with her Brother,
// filled with the primal froths of the holy Union,
// their little children, the eternal existance.


#endif
