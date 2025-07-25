function [optimizedList, resulGUI] = matRad_lexicographicOptimization(pln,dij,cst,priorityList)

    % Get the iteration scheme. This defines how the optimization sequence
    % will be carried out according to the list. In this case, this
    % iteration strategy simply goes through all the objectives.
    iteration = matRad_LexicographicIterationNominal(pln);

    % Initilaize the engine
    lexEngine = matRad_LexicographicEngine(iteration, priorityList, pln, cst, dij);

    % Optimize the list
    lexEngine.optimizeList();

    optimizedList = lexEngine.list;

    resulGUI = lexEngine.getResultsFromDir(dij, pln.propOpt.saveDir);
end