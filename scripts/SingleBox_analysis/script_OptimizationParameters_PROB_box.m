cst{1,6} = {};
cst{2,6} = {};
cst{3,6} = {};



cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(200,60));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));