// Scorer for TsEventsCounter
#include "TsEventsCounter.hh"

#include "G4SystemOfUnits.hh"

TsEventsCounter::TsEventsCounter(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{

	SetUnit("");

	//Load the energy binning specifications
	G4String name = GetFullParmName("Emax");
        fEmax = pM->GetDoubleParameter(name, "Energy");

	name = GetFullParmName("Emin");
        fEmin = pM->GetDoubleParameter(name, "Energy");

	name = GetFullParmName("nEBins");
        fEBins = pM->GetIntegerParameter(name);

	//Binning in Edep is turned off
	fEdepBins = 0;

        name = GetFullParmName("XBins");
        fPhantomBinningX = pM->GetIntegerParameter(name);

        name = GetFullParmName("ZBins");
        fPhantomBinningZ = pM->GetIntegerParameter(name);

        name = GetFullParmName("YBins");
        fPhantomBinningY = pM->GetIntegerParameter(name);


	fEBinning =  new G4double[fEBins + 1];
	fEdepBinning = new G4double[fEdepBins + 1];

	if (fPhantomBinningZ == (fEBins +1) && fPhantomBinningX == (fEdepBins +1)){

		GetBinning(fEmax,fEmin,fEBins,fEBinning);

	} else {
		printf("Warning, incorrect phantom binning. \n");
		printf("EBins = %i, Phantom Z Bins = %i \n", (fEBins +1) , fPhantomBinningZ);
		printf("EdepBins = %i, Phantom X Bins = %i \n", (fEdepBins +1) , fPhantomBinningX);
		for(int i=0;i<fEBins;i++){
			*(fEBinning+i) = 0;
		}
		for(int i=0;i<fEdepBins;i++){
			*(fEdepBinning+i) = 0;
		}
	}
}

TsEventsCounter::~TsEventsCounter()
{
	free (fEBinning);
	free (fEdepBinning);

}

G4bool TsEventsCounter::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

	ResolveSolid(aStep);

	G4double eDep = aStep->GetTotalEnergyDeposit();
	G4int particleAtomicMass = aStep->GetTrack()->GetDefinition()->GetAtomicMass();

	if (particleAtomicMass > 0 && eDep >0)
	{
		G4double totEPerNucleon = aStep->GetPreStepPoint()->GetKineticEnergy() / particleAtomicMass;

		//Get index of energy bin, add one so that Eindex = 0 is underflow
		G4int Eindex = interpolateIndex(totEPerNucleon, fEBinning, fEBins) + 1;

		//Get index of Edep
		G4int EdepIndex = 0;

		G4int DepthIndex = GetBin(GetIndex(aStep),1);

		long Idx = EdepIndex*(fEBins +1)*(fPhantomBinningY) + DepthIndex*(fEBins +1) + Eindex;

		G4double Weight = aStep->GetPreStepPoint()->GetWeight();
		AccumulateHit(aStep,Weight,Idx);
		//AccumulateHit(aStep,eDep,Idx);

	return true;
	}
return false;
}

G4int TsEventsCounter::interpolateIndex(G4double Value, G4double* Binning, G4int numberOfBins)
{
	G4int i = 0;
	while(*(Binning+i) < Value && i < numberOfBins)
		i++;
	i--;
	return i;
}

G4bool TsEventsCounter::GetBinning(G4double BinMax, G4double BinMin, G4int numOfBins, G4double* fBinning)
{
	G4double BinWidth = (BinMax - BinMin)/numOfBins;

	*(fBinning) = BinMin;

	G4int i = 0;

	while (*(fBinning+i) < BinMax && i < numOfBins){

		*(fBinning + i +1) = *(fBinning + i) + BinWidth;
		i++;
	};
	*(fBinning + i) = BinMax;

	return true;
}


