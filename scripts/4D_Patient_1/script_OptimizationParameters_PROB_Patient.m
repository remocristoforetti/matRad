for i=1:size(cst,1)
    cst{i,6} = {};
end



%cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,1));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,25));
cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,25));
cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,50));
