classdef Track_KC < Track
    properties
        Rc;
        Rp;
        Ek; %-> The total kinetic energy NOT Mev/u

        Z_ion;
        Beta_Ion;
        A;
        rho = 0.997; %g/cm^3
        Erest;
        AMU2MEV = 931.494;%938.272; %MeV/c^2 %Survival: 931.494027;
        Zeff;
        Kp;
        CONV = 0.1602;%Conversion factor for ((kev*cm^3)/(mu*g)) -> Gy*mu^2
    end

    methods
        function obj = Track_KC()
    
        end

        function dose = getDose(obj, r)
            %Compute track value as a function of radialDistance (in mu); 
           
            if (r < obj.Rc)
                dose = (1/(pi*obj.Rc^2))*(obj.CONV*(obj.LET/obj.rho) - 2*pi*obj.Kp*log(obj.Rp/obj.Rc));
            elseif (r < obj.Rp)
                dose = obj.Kp/(r^2);
            else
                dose = 0;
            end
        end

        function obj = updateTrack(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%% !!!!! %%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%% factor 2 in front added for testing
            %%%%%%%%%%%%%%%%%%%%%%%%% only!!!!
            obj.Rc = 0.0116*obj.Beta_Ion;                     %um
            obj.Rp = 0.0616*( (obj.Ek/obj.Erest)*obj.AMU2MEV )^1.7; %um
            
            %%%%%%%%%%%%%%%%%%%%%%%%% !!!!! %%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
 
            obj.Kp = (1.25*10^(-4))*(obj.Zeff/obj.Beta_Ion)^2;
            obj.radius = obj.Rp;
            obj.LET;
        end

        function Erest = get.Erest(obj)
            switch obj.Z_ion
                case 1 %p
                    A = 1;
                case 2 %He
                    A = 4;
                case 3  %Li
                    A = 7;
                case 4  %Be
                    A = 9;
                case 5  %B
                    A = 11;
                case 6  %C
                    A = 12;
            end
            obj.A = A;
            Erest = obj.A*obj.AMU2MEV;
        end

        function Beta_ion = get.Beta_Ion(obj)
            if ~isempty(obj.Ek) && ~isempty(obj.Erest)   
                Beta_ion = sqrt(1 - 1/((obj.Ek/obj.Erest +1)^2) );
            else
                Beta_ion = NaN;
            end
        end


        function Zeff = get.Zeff(obj)
            if ~isempty(obj.Z_ion) && ~isempty(obj.Beta_Ion)
                Zeff = obj.Z_ion*(1 - exp(-125*obj.Beta_Ion*(obj.Z_ion^(-2/3))));
            else
                Zeff = NaN;
            end
        end

        function obj = set.Ek(obj,value)
            %matRad_cfg = MatRad_Config();
            if value > 0
                obj.Ek = value;
            else
                %maRad_cfg.dispError('Negative Energy not supported. \n');
               disp('Negative Energy not supported. \n');
            end

            obj.updateTrack();
        end

        function obj = set.Z_ion(obj,value)
            %matRad_cfg = MatRad_Config();
            if (floor(value)==value) && (value < 7)
                obj.Z_ion = value;
            else
                %maRad_cfg.dispError('Z value not valid \n');
                disp('Z value not valid \n');
            end

            obj.updateTrack();
        end

        function plotTrack(obj,radii)
            if nargin<2
                radii = logspace(log10(obj.Rc)-2,log10(obj.Rp)); %um
            end

            LocalDose = arrayfun(@(x) obj.getDose(x), radii);

            figure;
            loglog(radii.*10^(-6), LocalDose, '.-');
        end

    end

end