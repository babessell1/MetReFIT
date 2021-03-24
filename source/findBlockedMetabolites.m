function [blockedMets,viableMets, modelT, demandRxns] = findBlockedMetabolites(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pull out intitial lists of blocked/unblocked metabolites by adding a 
%demand reaction and attempting to set as objective function for a model 
%where all media uptake reactions are allowed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
dbstop if error
%modelT = allowAllUptake(model);
modelT = model;
testMets = modelT.mets;
blockedMets = string([]);
viableMets = string([]);
for i = 1:length(testMets)
    tmet = testMets(i);
    disp(tmet)
    try
        modelT = addDemandReaction(modelT,tmet);
    catch 
        warning("Could not add demand reaction, would likely result in duplication");
    end
end

fprintf("Populating Demand Reaction List\n")
demandRxns = string([]);
for i = 1:length(modelT.rxns)
    rxn = modelT.rxns(i);
    if contains(rxn,"DM_")
        demandRxns = [demandRxns;rxn];
    end
end

blockedMets = findBlockedDemandReaction(modelT, demandRxns);

viableMets = setdiff(testMets, blockedMets);
blockedMets = sort(blockedMets);
viableMets = sort(viableMets);

toc
end

function blockedMetabolites = findBlockedDemandReaction(model, demandRxns)
% Determines those reactions which cannot carry any
% flux in the given simulation conditions.
%
% USAGE:
%
%    BlockedReaction = findBlockedReaction(model)
%
% INPUT:
%    model:               COBRA model structure
%
% OPTIONAL INPUT:
%    method:              'FVA' for flux variability analysis (default)
%                         'L2'  for 2-norm minimization
% OUTPUT:
%    blockedReactions:    List of blocked reactions
%
% .. Authors:
%       - Ines Thiele 02/09
%       - Srikiran C 07/14 - fixed error - assigning cells to blockedReactions which is a double
%       - Marouen BEN GUEBILA - used 2-norm min as a heuristic for non-sparsity
%       - Brandt Bessell 09/20 - made it exclusive to demand reactions

blockedReactions = string([]);
blockedMetabolites = string([]);
qptol = getCobraSolverParams('QP','intTol');
[m,n]=size(model.S);
% if (nargin < 2 || isequal(method, 'FVA'))
tol = 1e-10;
tic
demandRxns = cellstr(demandRxns);
[minMax(:, 1), minMax(:, 2)] = fluxVariability(model, 0, 'rxnNameList', demandRxns);
toc
cnt = 1;
for i = 1:length(minMax)
    if (minMax(i, 2) < tol && minMax(i, 2) > -tol && minMax(i, 1) < tol && minMax(i, 1) > -tol)
        blockedReactions(cnt) = demandRxns(i);
        cnt = cnt + 1;
    end
end
% else
%     model.c=zeros(n,1);
%     solution = solveCobraLPCPLEX(model, 0, 0, 0, [], 1e-6);
%     blockedReactions = model.rxns(abs(solution.full) < qptol)';
% end

for i = 1:length(blockedReactions)
    dm = string(blockedReactions(i));
    spl = strsplit(blockedReactions(i),"M_");
    blockedMetabolites = [blockedMetabolites;spl(2)]; 
end
end