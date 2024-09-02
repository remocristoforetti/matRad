function test_suite = test_FREDEngine

test_functions=localfunctions();

initTestSuite;

function test_getEngineFromPlnByName
    radModes = DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes;
    for i = 1:numel(radModes)
        plnDummy = struct('radiationMode',radModes{i},'machine','Generic','propDoseCalc',struct('engine','FRED'));
        engine = DoseEngines.matRad_ParticleFREDEngine.getEngineFromPln(plnDummy);
        assertTrue(isa(engine,'DoseEngines.matRad_ParticleFREDEngine'));
    end

function test_propertyAssignmentFromPln

    radModes = DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes;

    for i = 1:numel(radModes)
        pln = struct('radiationMode',radModes{i},'machine','Generic','propDoseCalc',struct('engine','FRED'));
        
        pln.propDoseCalc.HUclamping               = false;
        pln.propDoseCalc.HUtable                  = 'matRad_default_FRED';
        pln.propDoseCalc.exportCalculation        = true;
        pln.propDoseCalc.sourceModel              = 'emittance';
        pln.propDoseCalc.useWSL                   = true;
        pln.propDoseCalc.useGPU                   = false;
        pln.propDoseCalc.roomMaterial             = 'Vacuum';
        pln.propDoseCalc.printOutput              = false;
        pln.propDoseCalc.numHistoriesDirect       = 42;
        pln.propDoseCalc.numHistoriesPerBeamlet   = 42;
        
        engine = DoseEngines.matRad_ParticleFREDEngine.getEngineFromPln(pln);
        
        assertTrue(isa(engine,'DoseEngines.matRad_ParticleFREDEngine'));
        
        plnFields = fieldnames(pln.propDoseCalc);
        plnFields(strcmp([plnFields(:)], 'engine')) = [];
        
        for fieldIdx=1:numel(plnFields)
            assertTrue(isequal(engine.(plnFields{fieldIdx}), pln.propDoseCalc.(plnFields{fieldIdx})));
        end
    end

function test_writeFiles

        pln = struct('radiationMode','protons','machine','Generic','propDoseCalc',struct('engine','FRED'));



    