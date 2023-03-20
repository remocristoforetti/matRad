#ifndef TsEdEventsScorer_hh
#define TsEdEventsScorer_hh

#include "TsVBinnedScorer.hh"
#include "G4VProcess.hh"

class TsEdEventsScorer : public TsVBinnedScorer
{
public:
        TsEdEventsScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                        G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

        virtual ~TsEdEventsScorer();

protected:
	virtual G4int CombineSubScorers();
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) {return true;}
};
#endif
