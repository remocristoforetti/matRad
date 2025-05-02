classdef (Abstract) matRad_HPeffectVariance < matRad_ScalarQuantity
    % Just re-implementing some calculations here for simplicity and
    % clarity. This can be integrated later as a simpler class by defining
    % a subquantity for the mean effect squared
    properties (Abstract, Constant)
        
        vOmegaSubQuantity;
        effectSubQuantity;
    end

    methods
        function this = matRad_HPeffectVariance(cst)

            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});

        end

        function quantityOutput = computeQuantity(this, dij,struct,w)

            vOmegaSubQt = this.getSubQuantity(this.vOmegaSubQuantity);
            vMeanOmega  = vOmegaSubQt.getResult(dij,w);

            meanEffectSubQuantity = this.getSubQuantity(this.effectSubQuantity);
            
            % This will be a single distribution
            meanEffect            = meanEffectSubQuantity.getResult(dij,w);
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            currStructEffect = meanEffect{1}(currIdx);

            if ~iscolumn(currStructEffect)
                currStructEffect = currStructEffect';
            end

            quantityOutput = vMeanOmega{struct} - (1/N)*(currStructEffect'*currStructEffect);
            % quantityOutput = vMeanOmega{struct};
            % quantityOutput = -(1/N)*(currStructEffect'*currStructEffect);

            if quantityOutput<0
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Negative effect variance detected.');
            end
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            vOmegaSubQt    = this.getSubQuantity(this.vOmegaSubQuantity);
            vOmegaGradient = vOmegaSubQt.projectGradient(dij,struct,fGrad,w);

            meanEffectSubQuantity = this.getSubQuantity(this.effectSubQuantity);
            
            % This will be a single distribution
            meanEffect        = meanEffectSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            currStructEffect = zeros(size(meanEffect{1}));
            currStructEffect(currIdx) = meanEffect{1}(currIdx);


            fGradEffect = {2 .* (1/N) .*fGrad{struct}.* currStructEffect};

            meanEffectGrad = meanEffectSubQuantity.projectGradient(dij,1,fGradEffect,w);

            gradientOutput = vOmegaGradient - meanEffectGrad;
            % gradientOutput = vOmegaGradient;
            % gradientOutput = -meanEffectGrad;
        end
    end
end