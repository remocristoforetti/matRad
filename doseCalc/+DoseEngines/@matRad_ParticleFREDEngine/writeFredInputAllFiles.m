function writeFredInputAllFiles(this,stf)
    
    fred_cfg = MatRad_FREDConfig.instance();
    
    %write fred.inp file
    runFilename = fullfile(fred_cfg.MCrunFolder, fred_cfg.runInputFilename);
    this.writeRunFile(runFilename);

    %write region/region.inp file
    regionFilename = fullfile(fred_cfg.regionsFolder, fred_cfg.regionsFilename);
    this.writeRegionsFile(regionFilename, stf);

    %write plan file
    planFile = fullfile(fred_cfg.planFolder, fred_cfg.planFilename);
    this.writePlanFile(planFile,stf);

    %write fields file
    fieldsFile = fullfile(fred_cfg.planFolder,fred_cfg.fieldsFilename);
    this.writeFieldsFile(fieldsFile, stf);

    %write layers file
    layersFile = fullfile(fred_cfg.planFolder, fred_cfg.layersFilename);
    this.writeLayersFile(layersFile, stf);

    %write beamlets file
    beamletFile = fullfile(fred_cfg.planFolder, fred_cfg.beamletsFilename);
    this.writeBeamletsFile(beamletFile,stf);


    %write planDelivery file
    planFile = fullfile(fred_cfg.planFolder,fred_cfg.planDeliveryFilename);
    this.writeplanDeliveryFile(planFile, stf);
end
