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

/***************************************************************************
 * DESCRIPTION: Prisoner's Dilemma game engine by Marko Grönroos. 1997
 * Adopted from Mitchell's "An Introduction to Genetic Algorithms"
 ***************************************************************************/

#include <magic/mapplic.h>
#include <magic/mmath.h>
#include <nhp/genes.h>
#include <nhp/individual.h>
#include <nhp/simplepopula.h>
#include "prisoners.h"

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
 * Plays the IPD game with two participants.
 ******************************************************************************/
void PDGame::play (PDStrategy& plr0,  /**< Player object.           */
				   PDStrategy& plr1,  /**< Player object.           */
				   int         ngames /**< Number of games to play. */)
{
	mGameno = 0;
	mHistory.make (ngames, 2);

	double eq = (&plr0 == &plr1)? 2:1;

	for (int g=0; g<ngames; g++) {
		// Get decisions
		bool defect0 = plr0.getDecision (*this, 0);
		bool defect1 = plr1.getDecision (*this, 1);

		// TRACE2 ("0:%d    1:%d", int(defect0), int(defect1));

		// Calculate scores
		if ( defect0 &&  defect1)	{plr0.addScore (3/eq); plr1.addScore(3/eq);}
		if (!defect0 &&  defect1)	{plr0.addScore (5/eq); plr1.addScore(0/eq);}
		if ( defect0 && !defect1)	{plr0.addScore (0/eq); plr1.addScore(5/eq);}
		if (!defect0 && !defect1)	{plr0.addScore (1/eq); plr1.addScore(1/eq);}
		
		// Record history
		record (defect0, defect1);

		// Move on to next game
		nextGame ();
	}
}

/*******************************************************************************
 * Records result of a game iteration.
 ******************************************************************************/
void PDGame::record (bool defect0, /**< Did player 0 defect? */
					 bool defect1  /**< Did player 1 defect? */)
{
	mHistory.get (mGameno, 0) = defect0;
	mHistory.get (mGameno, 1) = defect1;
}

/*******************************************************************************
 * Returns the decision of player in a previous game.
 *
 * Throws an exception if trying to get a result before beginning.
 ******************************************************************************/
bool PDGame::getRecord (int gameno,
						int player) const
{
	ASSERT (gameno >= 0 && gameno < mGameno);
	ASSERT (player == 0 || player == 1);
	
	return mHistory.get (gameno, player);
}

/*******************************************************************************
 * \fn int PDGame::gameNo () const
 *
 * Returns the current successive game number.
 ******************************************************************************/

/*******************************************************************************
 * \fn void PDGame::nextGame ()
 *
 * Moves to next game.
 ******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           ----  ___    ----                                               //
//           |   ) |  \  (      |       ___   |   ___                        //
//           |---  |   |  ---  -+- |/\  ___| -+- /   )  ___  \   |           //
//           |     |   |     )  |  |   (   |  |  |---  (   \  \  |           //
//           |     |__/  ___/    \ |    \__|   \  \__   ---/   \_/           //
//                                                      __/   \_/            //
///////////////////////////////////////////////////////////////////////////////

void PDStrategy::make (const String& rules)
{
	ASSERTWITH (rules.length() == 70, "70 decisions required for a strategy");

	// First the 64-bit decision table of last three games
	for (int i=0; i<64; i++)
		mDecisionTable[i] = rules[i]=='1' || rules[i]=='D';

	// Then the hypothetical history
	for (int i=0; i<6; i++)
		mHypoHistory[i] = rules[i+64]=='1' || rules[i+64]=='D';
}

bool PDStrategy::getDecision (const PDGame& game,
							  int           player) const
{
	// Calculate history index
	int decisionIndex = 0;

	// Iterate from 3rd previous game to the previous game
	for (int round = game.gameNo()-3; round < game.gameNo(); round++) {
		int roundhistory = -1;

		// If record of the game is available
		if (round >= 0) {
			// Look at the history in the current game
			roundhistory =	game.getRecord (round, player) + 2*int(game.getRecord (round, 1-player));
		} else {
			// Otherwise the game has not been played, so we use hypothetical history
			roundhistory = mHypoHistory[(round+3)*2+player] + 2*mHypoHistory[(round+3)*2+1-player];
		}

		decisionIndex *= 4;
		decisionIndex += roundhistory;
	}

	// Return the decision according to the history index
	return mDecisionTable[decisionIndex];
}

TitForTatPDStrategy::TitForTatPDStrategy () {
	mName = "Tit4Tat";
	make ("0011001100110011001100110011001100110011001100110011001100110011000000");
}

RandomRulePDStrategy::RandomRulePDStrategy () {
	mName = "RandomRule";
	init ();
}

void RandomRulePDStrategy::init () {
	PDStrategy::init ();
	
	String str;
	str.reserve (70);
	for (int i=0; i<70; i++)
		str += (frnd()>0.5)? '1':'0';
	make (str);
}

bool RandomPDStrategy::getDecision (const PDGame& game, int player) const {
	return frnd()>0.5;
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ----      o                             ----                        //
//      |   )        ____        _    ___      |      ___    _    ___       //
//      |---  |/\ | (      __  |/ \  /   ) |/\ | --- /   ) |/ \  /   )      //
//      |     |   |  \__  /  \ |   | |---  |   |   \ |---  |   | |---       //
//      |     |   | ____) \__/ |   |  \__  |   |___/  \__  |   |  \__       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

impl_dynamic (PrisonerGene, {Gentainer});

/*******************************************************************************
* Adds private genes for the prisoner problem.
*
* The genome consists of an array of 70 binary genes.
*******************************************************************************/
void PrisonerGene::addPrivateGenes (Gentainer& g, const StringMap& params)
{
	for (int i=0; i<70; i++)
		add (new BinaryGene(format("D%d", i)));

	Gentainer::addPrivateGenes (g, params);
}

String PrisonerGene::translate () const
{
	String str;
	str.reserve (70);

	for (int i=0; i<70; i++)
		str += static_cast<const BinaryGene&>((*this)[i]).getvalue()? "1":"0";

	return str;
}

bool PrisonerGene::execute (const GeneticMsg& msg) const {
	msg.mrHost.set ("PDS", new String(translate ()));

	return false;
}

DataOStream& PrisonerGene::operator>> (DataOStream& out) const {
	out << translate ();
	return out; 
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ----      o                   ----   _   -----                    //
//        |   )        ____        _   |      / \  |       _                //
//        |---  |/\ | (      __  |/ \  | --- /   \ |---  |/ \  |   |        //
//        |     |   |  \__  /  \ |   | |   \ |---| |     |   |  \ /         //
//        |     |   | ____) \__/ |   | |___/ |   | |____ |   |   V          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

double evaluatePrisoner (
	const String& strDecs,
	bool          doPrint = false)
{
	int trials = 100,
		rounds = 100; // Games played at each trial

	// Make a league of contestants
	Array<PDStrategy> strategies;
	strategies.add (new PDStrategy(strDecs));      // The evolved strategy
	strategies.add (new TitForTatPDStrategy());
	strategies.add (new RandomRulePDStrategy());
	strategies.add (new RandomPDStrategy());
	strategies.add (new CooperatingPDStrategy());
	strategies.add (new DefectingPDStrategy());

	// Result table
	Vector scores (strategies.size());

	PDGame game;

	if (doPrint) {
		printf ("\n                ");
		for (int j=0; j<strategies.size(); j++)
			printf ("%1d(%8s)   ", j, (CONSTR) strategies[j].name());
		printf ("\n");
	}

	// For every strategy in league
	for (int i=0; i<strategies.size(); i++) {
		double total=0;
		if (doPrint)
			printf ("%1d(%10s): ", i, (CONSTR) strategies[i].name());

		// Compete strategy against every other strategy in the league
		for (int j=0; j<strategies.size(); j++) {

			// Play the game for a number of trials
			double s0=0, s1=0;
			for (int trial=0; trial<trials; trial++) {
				strategies[i].init ();
				strategies[j].init ();
				
				// Play a number of rounds
				game.play (strategies[i], strategies[j], rounds);
				
				// Record the result
				s0 += strategies[i].score()/rounds;
				s1 += strategies[j].score()/rounds;
			}
			s0 /= trials;
			s1 /= trials;

			if (doPrint)
				printf (" %2.3f/%2.3f  ", s0, s1);

			total += s0;
		}
		scores[i] = total/strategies.size();

		if (doPrint)
			printf ("  total=%2.3f\n", scores[i]);
	}

	return scores[0];
}

/*******************************************************************************
 * Adds environment specific features to genotype.
 ******************************************************************************/
void PrisonEAEnv::addFeaturesTo (
	Genome& genome /**< Genome to add the features to. */) const
{
	// Add the custom prisoner gene
	genome.add (new PrisonerGene ("PS"));
}

/*******************************************************************************
 * Evaluates individual (IPD strategy) by playing IPD against other strategies.
 ******************************************************************************/
double PrisonEAEnv::evaluateg (
	const Individual& ind /**< Individual to be evaluated in the environment. */)
{
	// Execute genotype to phenotype decoding for the "PS" (prisoner) gene.
	ind.execute (GeneticMsg ("PS", (Individual&) ind));
	
	// Read the phenotypic value from the individual
	const String& pds = static_cast<const String&> (ind["PDS"]);
	ASSERT (!isnull (pds));

	double score = evaluatePrisoner (pds, false);

	return score;
}

/*******************************************************************************
 * Writes evolution cycle report to output and log streams.
 *******************************************************************************/
void PrisonEAEnv::cycle_report (OStream& log, OStream& out)
{
	// Get decision table of the best individual in population
	String& pds = static_cast<String&> ((*mpBest)["PDS"]);
	ASSERT (!isnull (pds));

	// Evaluate the strategy
	double score = evaluatePrisoner (pds, true);

	printf ("score=%f\n", score);
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                        o                                  //
//                                   ___      _                              //
//                            |/|/|  ___| | |/ \                             //
//                            | | | (   | | |   |                            //
//                            | | |  \__| | |   |                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

Main ()
{
	sout.autoFlush ();
	readConfig ("prisoners.cfg");

	sout << toString(mParamMap) << "\n";

	try {
		sout << "Creating environment...\n";
		PrisonEAEnv prisonEnv;
		
		sout << "Creating population...\n";
		SimplePopulation pop (prisonEnv, mParamMap);

		sout << "Evolving...\n";
		pop.evolve (200, NULL);
	} catch (Exception& e) {
		sout << e.what();
	}
}
