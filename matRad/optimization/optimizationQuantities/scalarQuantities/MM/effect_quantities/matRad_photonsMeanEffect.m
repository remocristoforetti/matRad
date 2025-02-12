classdef matRad_photonsMeanEffect < matRad_MeanEffect

    properties (Constant)
        quantityName = 'photonsMeanEffect';
        requiredSubquantities = {'photonsMeanAlphaDose', 'photonsvMeanSqrtBeta'};

        alphaSubQuantity = 'photonsMeanAlphaDose';
        betaSubQuantity = 'photonsvMeanSqrtBeta';
    end

    methods
        function this = matRad_photonsMeanEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanEffect(supArg{:});
        end
    end
end