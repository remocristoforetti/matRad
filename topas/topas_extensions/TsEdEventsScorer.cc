// Scorer for TsEdEventsScorer

#include "TsEdEventsScorer.hh"
#include "TsParameterManager.hh"

TsEdEventsScorer::TsEdEventsScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                         G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("Gy");
	InstantiateSubScorer("TsEdEventsCounter", outFileName, "Numerator");
	InstantiateSubScorer("TsEventsCounter", outFileName, "Denominator");
}

TsEdEventsScorer::~TsEdEventsScorer() {;};

G4int TsEdEventsScorer::CombineSubScorers()
{
	TsVBinnedScorer* Numerator = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("Numerator"));
	TsVBinnedScorer* Denominator = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("Denominator"));

        for(unsigned index = 0; index<fFirstMomentMap.size(); index++)
        {
            fFirstMomentMap[index] = 0;
        };

	for(unsigned index = 0; index<fFirstMomentMap.size(); index++)
	{
		if(Denominator->fFirstMomentMap[index]>0)
			fFirstMomentMap[index] = Numerator->fFirstMomentMap[index]/Denominator->fFirstMomentMap[index];
	};

return 0;
}

