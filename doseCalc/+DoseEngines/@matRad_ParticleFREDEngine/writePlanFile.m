function writePlanFile(this, fName, stf)
    matRad_cfg = MatRad_Config.instance();

    fID = fopen(fName, 'w');
    

    try
        if this.calcDoseDirect
            simulatedPrimariesPerBixel = max([1, floor(this.numHistoriesDirect/sum([stf.totalNumOfBixels]))]);%this.numHistoriesPerBeamlet;%max([1, floor(this.numHistoriesDirect/stf.totalNumOfBixels)]);
        else
            simulatedPrimariesPerBixel = this.numHistoriesPerBeamlet;
        end
        
        fprintf(fID, 'nprim = %i\n', simulatedPrimariesPerBixel);


        layerCounter = 0;
        
        bixelCounter = 0;
        for i=1:numel(stf)
            for j=1:numel(stf(i).energies)
                fprintf(fID, '#Bixels Field%i, Layer%i\n', i-1,layerCounter+j-1);
                for k=1:stf(i).energyLayer(j).nBixels
                    currBixel.beamletID = num2str(bixelCounter+k-1);
                    currBixel.P         = arrayfun(@(idx) num2str(idx, '%2.3f'), [stf(i).energyLayer(j).rayPosX(k),stf(i).energyLayer(j).rayPosY(k),0], 'UniformOutput', false);
                    currBixel.v         = arrayfun(@(idx) num2str(idx, '%2.5f'), [stf(i).energyLayer(j).rayDivX(k),stf(i).energyLayer(j).rayDivY(k),1], 'UniformOutput', false);
                    currBixel.w         = num2str(stf(i).energyLayer(j).numOfPrimaries(k), '%2.7f');

                    printStructToDictionary(fID, currBixel, ['S', num2str(bixelCounter+k-1)],2);
                    
                end
                    
                    currLayer.Energy   = num2str(stf(i).energies(j));
                    currLayer.Espread  = num2str(stf(i).energySpread(j));
                    %currLayer.FWHM     = num2str(stf(i).FWHMs(j));
                    currLayer.emittanceX  = num2str(stf(i).emittanceX(j));
                    currLayer.twissAlphaX = num2str(stf(i).twissAlphaX(j));
                    currLayer.twissBetaX  = num2str(stf(i).twissBetaX(j));
                    
                    currLayer.beamlets = arrayfun(@(idx) ['S', num2str(idx)], bixelCounter:stf(i).energyLayer(j).nBixels+bixelCounter-1, 'UniformOutput', false);

                    
                    printStructToDictionary(fID, currLayer, ['L', num2str(layerCounter+j-1)],1);
                    fprintf(fID, '\n');
                    bixelCounter = bixelCounter + stf(i).energyLayer(j).nBixels;
            
            end
        
            
            fieldLimX = max(abs([stf(i).energyLayer.rayPosX])) + 10*max([stf(i).FWHMs]);
            fieldLimY = max(abs([stf(i).energyLayer.rayPosY])) + 10*max([stf(i).FWHMs]);

            currF.fieldNumber = i-1;
            currF.GA          = num2str(stf(i).gantryAngle);
            currF.CA          = num2str(stf(i).couchAngle);
            currF.ISO         = arrayfun(@num2str, stf(i).isoCenter, 'UniformOutput', false);
            currF.dim         = arrayfun(@num2str, [fieldLimX, fieldLimY, 0.1], 'UniformOutput', false);
            currF.Layers      = arrayfun(@(idx) ['L', num2str(idx)], layerCounter:numel(stf(i).energies)+layerCounter-1, 'UniformOutput', false);

            layerCounter = layerCounter + numel(stf(i).energies);

            
            printStructToDictionary(fID, currF, ['F', num2str(i-1)]);
            fprintf(fID, '\n');
        end

        %% Build plan
        plan.SAD    = stf(1).BAMStoIsoDist; 
        plan.Fields = arrayfun(@(i) ['F', num2str(i)], 0:numel(stf)-1, 'UniformOutput', false);
                
        printStructToDictionary(fID, plan, 'plan');
     
        
    catch
        matRad_cfg.dispError('Failed to write plan file');
    end

    fclose(fID);
end

function printStructToDictionary(fID, S, sName, indentTabs)

    if ~exist('indentTabs') || isempty(indentTabs)

        indentTabs = 0;
    end

    indentString = repmat('\t',1,indentTabs);
    
    fprintf(fID, indentString);
    
    fprintf(fID, 'def: %s = {', sName);

    sFields = fieldnames(S);

    for sFieldIdx =1:numel(sFields)

        currField = sFields{sFieldIdx};

        fprintf(fID, '''%s'': ',currField);

        if ~iscell(S.(currField))
            fprintf(fID, '%s', num2str(S.(currField)));

        else

            fprintf(fID, '[');
            for elementIdx=1:numel(S.(currField))-1
                fprintf(fID, '%s, ', S.(currField){elementIdx});
            end
            fprintf(fID, '%s]', S.(currField){end});

        end

        if sFieldIdx ~= numel(sFields)
            fprintf(fID, ', ');
        end

    end
    fprintf(fID, '}\n');

end