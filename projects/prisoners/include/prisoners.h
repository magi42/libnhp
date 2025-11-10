/***************************************************************************
 *   This file is part of the NeHeP library distribution.                  *
 *                                                                         *
 *   Copyright (C) 1998-2003 Marko Grönroos <magi@iki.fi>                  *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef __PRISONERS_H__
#define __PRISONERS_H__

#include <magic/mdatastream.h>
#include <nhp/gaenvrnmt.h>
#include <nhp/genetics.h>

/*******************************************************************************
 * Predeclarations
 ******************************************************************************/
class PDStrategy;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                    ----  ___    ----                                      //
//                    |   ) |  \  |      ___         ___                     //
//                    |---  |   | | ---  ___| |/|/| /   )                    //
//                    |     |   | |   \ (   | | | | |---                     //
//                    |     |__/  |___/  \__| | | |  \__                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Prisoner's dilemma game history recorder.
 *
 * Stores choises (cooperate/defect) made by each player in a single
 * iterative duel.
 ******************************************************************************/
class PDGame : public Object {
  public:
					PDGame		() {init();}

	void			init		() {mGameno=0;}
	void			play		(PDStrategy& plr0, PDStrategy& plr1, int ngames=1);
	bool			getRecord	(int gameno, int player) const;
	int				gameNo		() const {return mGameno;}
	
  protected:
	void			nextGame	() {mGameno++;}
	void			record		(bool defect0, bool defect1);

  private:
	PackTable<int>	mHistory;
	int				mGameno;
};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           ----  ___    ----                                               //
//           |   ) |  \  (      |       ___   |   ___                        //
//           |---  |   |  ---  -+- |/\  ___| -+- /   )  ___  \   |           //
//           |     |   |     )  |  |   (   |  |  |---  (   \  \  |           //
//           |     |__/  ___/    \ |    \__|   \  \__   ---/   \_/           //
//                                                      __/   \_/            //
///////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Iterated Prisoner's Dilemma Strategy.
 *
 * Every player plays PD according to a strategy defined in a class
 * inherited from this baseclass.
 *
 * The strategy is defined as a decision table given in the
 * constructor for this class. Alternatively, an inheritor can
 * reimplement @ref getDecision() to provide a custom decision scheme.
 *
 * This baseclass handles also scoring.
 ******************************************************************************/
class PDStrategy : public Object {
  public:
						PDStrategy	() : mName ("Any") {
							mDecisionTable.make (70);
							mHypoHistory.make (6);
							mScore = 0;
						}
						PDStrategy	(const String& rules) : mName ("Any") {
							mDecisionTable.make (70);
							mHypoHistory.make (6);
							mScore = 0;
							make (rules);
						}

	void				make		(const String& rules);

	// Initializes the strategy (sets score to 0)
	virtual void		init		() {mScore=0;}
	
	// Returns a decision for the given game situation
	virtual bool		getDecision	(const PDGame& game, int playerRole) const;

	void				addScore	(double x) {mScore += x;}

	double				score		() const {return mScore;}
	const String&		name		() const {return mName;}
	
	virtual PDStrategy*	clone		() {return new PDStrategy ();}

  protected:
	String			mName;

  private:
	PackArray<int>	mHypoHistory;    /**< Hypothetical history. */
	PackArray<int>  mDecisionTable;
	double			mScore;
};

// Do what the other did in the last game. Presume cooperation
class TitForTatPDStrategy : public PDStrategy {
  public:
					TitForTatPDStrategy	();
	PDStrategy*		clone		() {return new TitForTatPDStrategy ();}
};

// Play with a random rule set
class RandomRulePDStrategy : public PDStrategy {
  public:
					RandomRulePDStrategy	();
	virtual void	init		();
	PDStrategy*		clone		() {return new RandomRulePDStrategy ();}
};

// Play totally randomly
class RandomPDStrategy : public PDStrategy {
  public:
					RandomPDStrategy	() {mName = "XRandom";}
	virtual bool	getDecision	(const PDGame& game, int player) const;
	PDStrategy*		clone		() {return new RandomPDStrategy ();}
};

// Play totally randomly
class CooperatingPDStrategy : public PDStrategy {
  public:
					CooperatingPDStrategy	() {mName = "Cooping";}
	virtual bool	getDecision	(const PDGame& game, int player) const {return false;} 
	PDStrategy*		clone		() {return new CooperatingPDStrategy ();}
};

// Play totally randomly
class DefectingPDStrategy : public PDStrategy {
  public:
					DefectingPDStrategy	() {mName ="Defecting";}
	virtual bool	getDecision	(const PDGame& game, int player) const {return true;}
	PDStrategy*		clone		() {return new DefectingPDStrategy ();}
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ----      o                             ----                        //
//      |   )        ____        _    ___      |      ___    _    ___       //
//      |---  |/\ | (      __  |/ \  /   ) |/\ | --- /   ) |/ \  /   )      //
//      |     |   |  \__  /  \ |   | |---  |   |   \ |---  |   | |---       //
//      |     |   | ____) \__/ |   |  \__  |   |___/  \__  |   |  \__       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class PrisonerGene : public Gentainer {
	decl_dynamic (PrisonerGene);
  public:
						PrisonerGene	(const GeneticID& name = NULL) : Gentainer (name) {;}
						PrisonerGene	(const PrisonerGene& orig) : Gentainer (orig) {;}

	String				translate		() const;
	
	virtual void		addPrivateGenes	(Gentainer& g, const StringMap& params);
	virtual bool		execute			(const GeneticMsg& msg) const;
	virtual Genstruct*	replicate		() const {return new PrisonerGene (*this);}
	DataOStream&		operator>>		(DataOStream& out) const; 
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ----      o                   ----   _   -----                    //
//        |   )        ____        _   |      / \  |       _                //
//        |---  |/\ | (      __  |/ \  | --- /   \ |---  |/ \  |   |        //
//        |     |   |  \__  /  \ |   | |   \ |---| |     |   |  \ /         //
//        |     |   | ____) \__/ |   | |___/ |   | |____ |   |   V          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class PrisonEAEnv : public EAEnvironment {
  public:

	virtual void	addFeaturesTo	(Genome& genome) const;
	virtual double	evaluateg		(const Individual& ind);
	virtual void	cycle_report	(OStream& log, OStream& out);
};

#endif
