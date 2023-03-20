// Scorer for TsDoseEventsCounter
#include "TsDoseEventsCounter.hh"

#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"

TsDoseEventsCounter::TsDoseEventsCounter(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fEmCalculator()
{
        //SetUnit("Gy");

        fWater = GetMaterial("G4_WATER");

        fProtonSubstituteForNeutrals = G4Proton::ProtonDefinition();

        if (fPm->ParameterExists(GetFullParmName("SubstituteEnergyForNeutralScaling")))
                fSubstituteForNeutralsEnergy =  fPm->GetDoubleParameter(GetFullParmName("SubstituteEnergyForNeutralScaling"), "Energy");
        else
                fSubstituteForNeutralsEnergy =  100*MeV;
//}

//{
	//The subscorer inherits the binng of the scorer, somehow this still works and load the correct binning. This part with the load for the parameters should go into the sccorer that calls this subscorer.
	SetUnit("Gy");

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

TsDoseEventsCounter::~TsDoseEventsCounter()
{
	free (fEBinning);
	free (fEdepBinning);

}

G4bool TsDoseEventsCounter::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

	ResolveSolid(aStep);

	G4double eDep = aStep->GetTotalEnergyDeposit();
	G4int particleAtomicMass = aStep->GetTrack()->GetDefinition()->GetAtomicMass();

	if (particleAtomicMass > 0 && eDep >0)
	{
		G4double totEPerNucleon = aStep->GetPreStepPoint()->GetKineticEnergy() / particleAtomicMass;

		// This dose is already mult. by weight, is dose per single event
		G4double dose = getDose(aStep, eDep);

		//Get index of energy bin, add one so that Eindex = 0 is underflow
		G4int Eindex = interpolateIndex(totEPerNucleon, fEBinning, fEBins) + 1;

		//Get index of Edep
		G4int EdepIndex = 0;

		G4int DepthIndex = GetBin(GetIndex(aStep),1);

		long Idx = EdepIndex*(fEBins +1)*(fPhantomBinningY) + DepthIndex*(fEBins +1) + Eindex;

		//G4double Weight = aStep->GetPreStepPoint()->GetWeight();

		AccumulateHit(aStep,dose,Idx);
		//AccumulateHit(aStep,eDep,Idx);

	return true;
	}
return false;
}

G4int TsDoseEventsCounter::interpolateIndex(G4double Value, G4double* Binning, G4int numberOfBins)
{
	G4int i = 0;
	while(*(Binning+i) < Value && i < numberOfBins)
		i++;
	i--;
	return i;
}

G4bool TsDoseEventsCounter::GetBinning(G4double BinMax, G4double BinMin, G4int numOfBins, G4double* fBinning)
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


G4double TsDoseEventsCounter::getDose(G4Step* aStep, G4double eDep)
{
	if ( eDep > 0. ) {
		G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();

		ResolveSolid(aStep);

		G4double dose = eDep / ( density * fSolid->GetCubicVolume() );

		dose *= aStep->GetPreStepPoint()->GetWeight();

		G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();

		if ( particle->GetPDGCharge() != 0 ) {
			// convert to Dose to Water for charged particles from EM tables:
			G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
			G4double materialStoppingPower = fEmCalculator.ComputeTotalDEDX(energy, particle,aStep->GetPreStepPoint()->GetMaterial());

			G4double waterStoppingPower = fEmCalculator.ComputeTotalDEDX(energy, particle, fWater);
			dose *= ( density / fWater->GetDensity() ) *  ( waterStoppingPower / materialStoppingPower );
		} else {
			// convert to Dose to Water for neutral particles from EM tables, assuming scaling factor from protons (dE/dx not defined for neutrals):
			G4double materialStoppingPower = fEmCalculator.ComputeTotalDEDX(fSubstituteForNeutralsEnergy, fProtonSubstituteForNeutrals,
																			aStep->GetPreStepPoint()->GetMaterial());
			G4double waterStoppingPower = fEmCalculator.ComputeTotalDEDX(fSubstituteForNeutralsEnergy, fProtonSubstituteForNeutrals, fWater);
			dose *= ( density / fWater->GetDensity() ) *  ( waterStoppingPower / materialStoppingPower );
		}
	return dose;
	};
return 0;
}






