// Statistical module.
// Collects statistics about mutations rates.

// The globality is "a bit ugly"

#ifndef __MUTRECORD_H__
#define __MUTRECORD_H__

/** A logging class that stores some statistical data during evolution
 *  runs.
 **/
class MutabilityRecord : public Object {
  public:

	/** Is the recording on or off.
	 **/
	static bool		record;

	/** Adds one statistical instance to @ref FloatGene's mutation
	 *  rate record.
	 **/
	static void		addFloatMutability	(double x) {
		smFloatSum+=x;
		smFloatSamples++;
		if (x<smFloatMin) smFloatMin = x;
		if (x>smFloatMax) smFloatMax = x;
	}
	/** Adds one statistical instance to @ref FloatGene's mutation
	 *  variance record.
	 **/
	static void		addFloatVariance	(double x) {
		smFloatVarSum+=x;
		smFloatVarSamples++;
		if (x<smFloatVarMin) smFloatVarMin = x;
		if (x>smFloatVarMax) smFloatVarMax = x;
	}
	/** Adds one statistical instance to @ref BooleanGene's mutation
	 *  rate record.
	 **/
	static void		addBoolMutability	(double x) {
		smBoolSum+=x;
		smBoolSamples++;
		if (x<smBoolMin) smBoolMin = x;
		if (x>smBoolMax) smBoolMax = x;
	}
	/** Resets the stats.
	 **/
	static void		reset				() {
		smFloatVarMin=smFloatMin=smBoolMin=666;
		smFloatVarSum=smFloatSum=smBoolSum=0;
		smFloatVarMax=smFloatMax=smBoolMax=-1;
		smFloatVarSamples=smFloatSamples=smBoolSamples=0;
	}

	static double	floatMin			() {return smFloatMin;}
	static double	floatAvg			() {return smFloatSum/smFloatSamples;}
	static double	floatMax			() {return smFloatMax;}
	static double	floatVarMin			() {return smFloatVarMin;}
	static double	floatVarAvg			() {return smFloatVarSum/smFloatVarSamples;}
	static double	floatVarMax			() {return smFloatVarMax;}
	static double	boolMin				() {return smBoolMin;}
	static double	boolAvg				() {return smBoolSum/smBoolSamples;}
	static double	boolMax				() {return smBoolMax;}
	
  protected:
	static int		smFloatVarSamples, smFloatSamples, smBoolSamples;
	static double	smBoolSum, smBoolMin, smBoolMax;
	static double	smFloatSum, smFloatMin, smFloatMax;
	static double	smFloatVarSum, smFloatVarMin, smFloatVarMax;
};

#endif
