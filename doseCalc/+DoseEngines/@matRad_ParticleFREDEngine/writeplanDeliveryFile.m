function writeplanDeliveryFile(this, fName, stf)
    matRad_cfg = MatRad_Config.instance();
    fred_cfg = MatRad_FREDConfig.instance();

    %Try to read and copy template for pln delivery
    try
        planDeliveryTemplateFilename = fullfile(pwd, fred_cfg.planDeliveryTemplate);
        templateFile = fileread(planDeliveryTemplateFilename);

        fID = fopen(fName, 'w');
        fprintf(fID, templateFile);
        fclose(fID);

    catch
        matRad_cfg.dispWarning('Unable to read plan delivery template, writing one form scratch. This might not be uptated!')
    

        fID = fopen(fName, 'w');
        try
            fprintf(fID, 'include: inp/plan/fields.inp\n');
            fprintf(fID, 'def: ipb = 0\n');
            
            
            %deactivate fileds
            fprintf(fID, 'for(fieldIdx in range(nFields))<\n');
                fprintf(fID, '\tdeactivate: field_$fieldIdx\n');
            fprintf(fID, 'for>\n\n');

                %loop over fields
                fprintf(fID, 'for(fieldIdx in range(nFields))<\n');
                    
                
                %activate current filed
                    fprintf(fID, '\tactivate: field_$fieldIdx\n');
                    fprintf(fID, '\tdef: GA = valueArr(PatientGA,fieldIdx)\n');
                    fprintf(fID, '\tdef: CA = valueArr(PatientCA,fieldIdx)\n');
                    fprintf(fID, '\tdef: ISO = getDictionaryKeyValue(PatientISO,"ISOField_"+str(fieldIdx))\n');

                    fprintf(fID, '\ttransform: Phantom move_to ${valueArr(ISO,0)} ${valueArr(ISO,1)} ${valueArr(ISO,2)} Room\n');
                    fprintf(fID, '\ttransform: Phantom rotate y ${-1*CA} Room\n');
                    fprintf(fID, '\ttransform: Phantom rotate z ${-1*GA} Room\n\n');
                    
                    %define pbmaster
                    fprintf(fID, '\tpbmaster: $fieldIdx; Xsec = gauss; columns = [P.x, P.y, N, FWHM, T, v.x, v.y, v.z]\n\n');
                    fprintf(fID, '\tdef: layersInField = getDictionaryKeyValue(layers,"Field_"+str(fieldIdx))\n\n');
                    fprintf(fID, '\tfor(layer in layersInField)<\n\n');
                    fprintf(fID, '\t\tdef: currEnergy = getDictionaryKeyValue(layer,"Energy")\n');
                    fprintf(fID, '\t\tdef: currFWHM   = getDictionaryKeyValue(layer,"FWHM")\n\n');
                    fprintf(fID, '\t\tfor(pbIdx in range(layer["nSpots"]))<\n');
                        fprintf(fID, '\t\t\tdef: x = valueArr(getDictionaryKeyValue(layer,"x"),pbIdx)\n');
                        fprintf(fID, '\t\t\tdef: y = valueArr(getDictionaryKeyValue(layer,"y"),pbIdx)\n');
                        fprintf(fID, '\t\t\tdef: w = %e*valueArr(getDictionaryKeyValue(layer,"w"),pbIdx)\n', this.conversionFactor);
                        fprintf(fID, '\t\t\tdef: vx = valueArr(getDictionaryKeyValue(layer,"divX"),pbIdx)\n');
                        fprintf(fID, '\t\t\tdef: vy = valueArr(getDictionaryKeyValue(layer,"divY"),pbIdx)\n\n');

                        fprintf(fID, '\t\t\tpb: $ipb $fieldIdx $x $y $w $currFWHM $currEnergy $vx $vy 1\n\n');

                        fprintf(fID, '\t\t\tdef: ipb = ipb +1;\n');
                    fprintf(fID, '\t\tfor>\n');
                fprintf(fID, '\tfor>\n');
                fprintf(fID, '\tdeliver: field_$fieldIdx\n');
                fprintf(fID, '\tdeactivate: field_$fieldIdx\n');
                %fprintf(fID, '\ttransform: Phantom rotate z ${-1*GA} self\n');
                fprintf(fID, '\ttransform: Phantom rotate z ${GA} Room\n');
                fprintf(fID, '\ttransform: Phantom rotate y ${CA} Room\n');
                fprintf(fID, '\ttransform: Phantom move_to 0 0 0 Room\n\n');

                fprintf(fID, 'for>\n\n');
            
        catch
            matRad_cfg.dispError('Failed to write planDelivery file');
        end

        fclose(fID);
    end

end

