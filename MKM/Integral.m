classdef Integral < handle
    properties
        r1 = 0;
        r2 = 1e-4;
        log10r2;
        y= 0;
        z1 = 0;
        z2 = 0;
        stepWeightedIntegral = 0.01;
        corFac = 2;
    end


    methods
        function obj = Integral()
            %Constructor;
        end


        
        function obj.computeZD(obj, Track, domain, Resolution)

                r_max = Track.radius;
                rNuc = domain.rNuc;
                nSteps = ceil((r_max + rNuc)/Resolution);
                y = linspace(0,r_max + rNuc,nSteps);

                d2 = arrayfun(@(x) obj.computeDose(x, rNuc, Track), y);
                dr = y(1:end-1) - y(2:end);
                Z1 = sum((d2(1:end-1).*y(1:end-1) + d2(2:end).*y(2:end)).*(y(1:end-1) - y(2:end)));
                Z1D = sum(((d2(1:end-1).^2).*y(1:end-1) + (d2(2:end).^2).*y(2:end)).*(y(1:end-1) - y(2:end)));
                
                obj.gamma = Z1D/Z1;
                obj.E = Track.E;
                
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
                    area1 = pi*r_intresection^2;
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
                            area3 = pi*rNuc^2 - (obj.corFact*theta1*rMax^2 + obj.corFact*theta2*rNuc^2 - rMax*y*sin(theta1));
                        else
                            area3 =pi*(rNuc^2 - r_intersection^2);
                        end
                    end
                elseif (y <= (rNuc + Track.radius)) %if still intesection, but not getting track center
                    rMin = y - rNuc;
                    [d2,area2] = Track.WeightedIntegral(rMin,rMax,y,obj.stepWeightedIntegral);

                    if (rMax == Track.radius)
                        theta1 = acos( (y^2 + rMax^2 - rNuc^2)/(2*y*rMax) );
                        theta2 = acos( (y^2 - rMax^2 + rNuc^2)/(2*y*rNuc) );
                        area3 = pi*rNuc^2 - (obj.corFact*theta1*rMax^2 + obj.corFact*theta2*rNuc^2 - rMax*y*sin(theta1));
                    end
                else


                end
                dose_y = (d1 + d2)/(area1 + area2 +area3);

        end
    end
end