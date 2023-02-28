classdef IntegralDose < handle
    properties
        gamma;
        E;
        stepWeightedIntegral = 0.01;

        corFact = 1;
        Resolution = 0.01;
        Zsat;
    end


    methods
        function obj = IntegralDose()
            %Constructor;
        end


        
        function computeZD(obj, Track, domain)

                r_max = Track.radius;
                rNuc = domain.rNuc;
%                 nSteps = ceil((r_max + rNuc)/obj.Resolution);
%                 y = [0, logspace(log10(obj.Resolution),log10(r_max + rNuc),nSteps)];
                %y = logspace(-3,3,1000);
                %Step = 10^(-3);
                Step = obj.Resolution;
                y(1) = 0;
                y(2) = 10^(-4);
                k=2;
                while(y(k)<(r_max + rNuc))
                    y(k+1) = 10^(log10(y(k) + Step));
                    k = k+1;
                end
                d2 = arrayfun(@(x) obj.computeDose(x, rNuc, Track), y);
                
                dr = -(y(1:end-1) - y(2:end));
                Z1 = sum((d2(1:end-1).*y(1:end-1) + d2(2:end).*y(2:end)).*dr);
                %Z1D = sum(((d2(1:end-1).^2).*y(1:end-1) + (d2(2:end).^2).*y(2:end)).*(y(1:end-1) - y(2:end)), 'omitNaN');

                RadiiRatioSquared = (domain.RN/rNuc)^2;
                Z0 = RadiiRatioSquared/(sqrt(domain.Beta_Tissue*(1+ RadiiRatioSquared)));
                Zsat = ((Z0^2)./d2).*(1-exp(-((d2./Z0).^2)));
                ZatD2 = Zsat.*d2;


                Z1D = sum((ZatD2(1:end-1).*y(1:end-1) + ZatD2(2:end).*y(2:end)).*(dr), 'omitNaN');
                
                obj.gamma = Z1D/Z1;
                obj.E = Track.Ek;
                %obj.Zsat = Zsat;
                %-> define y (array) -> ImpactParameter
                %d2(y) = arrayfun(@(y) computeDose(y))
                %z1(d2)
                %zD1(d2)
                %Zd
        end
        
        function dose_y = computeDose(obj, y, rNuc, Track)
                area1 = 0;
                area2 = 0;
                area3 = 0;
                d1 = 0;
                d2 = 0;
                rMax = min([Track.radius,y + rNuc]);
                
                if (y < rNuc)
                    if (y + Track.radius < rNuc)
                        r_intersection = Track.radius;
                    else
                        r_intersection = rNuc - y;
                    end

                    %Compute A1 and D1
                    area1 = pi*r_intersection^2;
                    d1 = Track.radialIntegral(0,r_intersection)*area1;
                    
                    %Compute D2 and A2    
                    if (rMax > r_intersection)
                       [td2, area2] = Track.WeightedIntegral(r_intersection, rMax, y, obj.stepWeightedIntegral,rNuc);
                       d2 = td2*area2;
                    end

                    %Compute A3 (no dose)
                    if (rMax == Track.radius)
                        if (Track.radius > rNuc - y)
                            theta1 = acos( (y^2 + rMax^2 - rNuc^2)/(2*y*rMax) );
                            theta2 = acos( (y^2 - rMax^2 + rNuc^2)/(2*y*rNuc) );
                            area3 = pi*rNuc^2 - (obj.corFact*theta1*(rMax^2) + obj.corFact*theta2*(rNuc^2) - rMax*y*sin(theta1));
                        else
                            area3 =pi*(rNuc^2 - r_intersection^2);
                        end
                    end
                elseif (y <= (rNuc + Track.radius)) %if still intesection, but not getting track center
                    rMin = y - rNuc;
                    [d2t,area2] = Track.WeightedIntegral(rMin,rMax,y,obj.stepWeightedIntegral, rNuc);
                    d2 = d2t*area2;
                    if (rMax == Track.radius)
                        arg = (y^2 + rMax^2 - rNuc^2)/(2*y*rMax);
                        if abs(arg) >1
                            arg = 1;
                        end
                        theta1 = acos( arg );
                        
                        arg = (y^2 - rMax^2 + rNuc^2)/(2*y*rNuc);
                        if abs(arg) >1
                            arg = 1;
                        end
                        theta2 = acos( arg );
                        area3 = pi*rNuc^2 - (obj.corFact*theta1*rMax^2 + obj.corFact*theta2*rNuc^2 - rMax*y*sin(theta1));
                    end
                else

                d1 = 0;
                d2 = 0;
                end
                if sum(area1 + area2 + area3) >0
                    dose_y = (d1 + d2)/(area1 + area2 +area3);
                else
                    dose_y = 0;

                end

        end

        function plotzD(obj,track,domain)
            if nargin <3
                track = Track_KC();
                track.Z_ion = 6;
                track.Ek = 50*track.A;
                domain = Domain();
                domain.rNuc = 0.32;

            end

            rNuc = domain.rNuc;
            r_max = track.radius;

            if isempty(obj.Resolution)
                y = logspace(-3,3,1000);
            else
                Step = obj.Resolution;
                y(1) = 0;
                y(2) = 10^(-4);
                k=2;
                while(y(k)<(r_max + rNuc))
                    y(k+1) = 10^(log10(y(k) + Step));
                    k = k+1;
                end
            end
            
            d2 = arrayfun(@(x) obj.computeDose(x, rNuc, track), y);

            %figure;
            
            loglog(y*10^(-6), d2, '.-');
            grid on;
            xlabel('Impact Parameter [m]');
            ylabel('Specific Energy [Gy]');
            title('Carbon ion, E = ', (track.Ek/track.A));
            xlim([10^(-9), 10^(-4)]);
            hold on;
        end
    end
end