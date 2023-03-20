#ifndef TsEdEventsCounter_hh
#define TsEdEventsCounter_hh

#include "TsVBinnedScorer.hh"
#include "G4VProcess.hh"

class TsEdEventsCounter : public TsVBinnedScorer
{
public:
	TsEdEventsCounter(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
			G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~TsEdEventsCounter();

	G4bool GetBinning(G4double BinMax, G4double BinMin, G4int numOfBins, G4double* fBinning);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    	G4int interpolateE(G4double EnergyPerNucleon);
	G4int interpolateIndex(G4double Value, G4double* Binning, G4int numOfBins);

protected:

    G4double fEmax;
    G4double fEmin;
    G4int fEBins;

    G4double fEdepMax;
    G4double fEdepMin;
    G4int fEdepBins;

    G4double* fEBinning;
    G4double* fEdepBinning;

    G4int fPhantomBinningX;
    G4int fPhantomBinningZ;
    G4int fPhantomBinningY;
};
#endif
