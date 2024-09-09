classdef matRad_PhotonStfGenerator < matRad_ExternalStfGenerator

    properties (Constant)
        name = 'photonStfGen';
        shortName = 'photonStfGen';
        possibleRadiationModes = {'photons'};
    end 

    
    
    methods 
        function this = matRad_PhotonStfGenerator(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_ExternalStfGenerator(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'photons';
            end
        end            
    end

    methods (Access = protected)        
        function pbMargin = getPbMargin(this)
            pbMargin = this.bixelWidth;
        end
        function rayPos = initializeRayTargetPosition(this,rayPos,rotMat_vectors_T,SAD) 

            rayPos = this.initializeRayTargetPosition@matRad_ExternalStfGenerator(rayPos,rotMat_vectors_T,SAD);

            %photon ray-target position
            for j = 1:rayPos.numOfRays
                        rayPos.ray(j).beamletCornersAtIso = [this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2]]*rotMat_vectors_T;
                        rayPos.ray(j).rayCorners_SCD = (repmat([0, this.machine.meta.SCD - SAD, 0],4,1)+ (this.machine.meta.SCD/SAD) * ...
                            [this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2]])*rotMat_vectors_T;
            end
        end
        function rays = initializeEnergy(this,rays,ct)
            for j = rays.numOfRays:-1:1

                rays.ray(j).energy = this.machine.data.energy;
            end
        end
    end
end
