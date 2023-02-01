classdef Track < handle
    properties
       radius;
       A_1 = 4.382*10^(-25); %(8.4018*10^(-54))/(4*pi*0.51); n*e^4/4*pi*mec^2*e0^2; %https://www.nature.com/articles/s41598-017-10554-0
       A_2 = 0.511; %mec^2

       A_3 = 3.893*10^(22); %Conversion factor (6.24*10^9)^2 * 10^3 * 10^(-6);
       LET; %Kev/mu
       I = 75*10^(-6); %MeV

%        n = 3.343*10^(29);
%        e = 1.602*10^(-19);
%        e0 = 8.854*10^(-12);
    end


    methods
        function obj = Track();

        end

        function dose = radialIntegral(obj,rMin,rMax)
                %Computes radial integral of Track structure (as a function of radial distance)
               % obj.updateTrack();
                rEnd = min([rMax, obj.Rp]);
                if rMin < obj.Rc
                    if rMax < obj.Rc
                        dose = obj.getDose(rEnd)*(rMax^2 - rMin^2)/2; %-> Constant central value from central interal. Integral over angle gives always 2pi and then get simplified because of normalization. So it is omitted.
                    elseif rMax <= obj.Rp
                        dose = obj.getDose(rMin)*(obj.Rc^2 - rMin^2)/2 + 2*obj.Kp*log(rEnd/obj.Rc);
                    end
                elseif rMin < obj.Rp
                    dose = 2*obj.Kp*log(rEnd/rMin);

                else
                    dose = 0;
                end
                dose = dose/(rMax^2 -rMin^2);
        end

        function [dose, area] = WeightedIntegral(obj,rMin,rMax,y,step, rNuc)
                 nSteps = ceil((rMax-rMin)/step);
                %%%%%%%%%%%%%% !!!!!%%%%%%%%%%%%%% Changed nSteps, was 10
                 if nSteps <10
                    %step = (log10(rMax) -log10(rMin))/10;
                    nSteps = 10;
                end

                %Steps = [rMin,10.^(log10(rMin) + step*[1:nSteps])];
                if nSteps > 0
                    Steps = logspace(log10(rMin),log10(rMax),nSteps);
                    W = arrayfun(@(r) obj.getArcWeight(r,y,rNuc), Steps);
                    F = arrayfun(@(r) obj.getDose(r),Steps);
                    F = F.*Steps.*W;
                    dr = Steps(2:end) - Steps(1:end-1);
                    InteGral = sum(((F(2:end) + F(1:end-1))/2).*dr);
                    area = sum(((W(1:end-1).*Steps(1:end-1) + W(2:end).*Steps(2:end))/2).*dr);
                    if area>0
                        dose = InteGral/area;
                    else
                        dose = 0;
                    end
                else
                    dose = 0;
                    area = 0;
                end
        end


        function w = getArcWeight(obj,r,y,rNuc)
              if (y < rNuc) 
                    if (r <= rNuc - y)
                        w = 2*pi;
                    elseif (r < y + rNuc)
                        arg = y/(2*r) + r/(2*y) - rNuc*rNuc / (2*y*r);
                        if abs(arg)>1
                            arg = 1;
                        end
                        w = 2 * acos( arg );
                    else
                        w=0;
                    end
              else
                    if (r <= y - rNuc)
                        w =  0;
                    elseif (r <= y + rNuc)
                        arg = y/(2*r) + r/(2*y) - rNuc*rNuc / (2*y*r);
                        if (arg > 1.)
                            arg = 1;
                        end
                        w= 2*acos( arg );
                    else
                        w=0;

                    end
              end
%             arg = y/(2*r) + r/(2*y) - rNuc*rNuc/(2*y*r);
%             if abs(arg) > 1%(r <= y + rNuc) && (abs(arg) > 1)
%                 arg = 1;
%             end
%             w = 2*acos(arg);
        end

        function LET = get.LET(obj)
            LOG = log((2*obj.A_2/obj.I)*((obj.Beta_Ion^2)/(1-obj.Beta_Ion^2)));
            LET = obj.A_1*((obj.Z_ion/obj.Beta_Ion)^2)*(LOG - obj.Beta_Ion^2)*obj.A_3;
        end
    end
end
