// Scorer for TsDoseEventsScorer

#include "TsDoseEventsScorer.hh"
#include "TsParameterManager.hh"

TsDoseEventsScorer::TsDoseEventsScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                         G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("Gy");
	InstantiateSubScorer("TsDoseEventsCounter", outFileName, "Numerator");
	InstantiateSubScorer("TsEventsCounter", outFileName, "Denominator");
}

TsDoseEventsScorer::~TsDoseEventsScorer() {;};

G4int TsDoseEventsScorer::CombineSubScorers()
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

