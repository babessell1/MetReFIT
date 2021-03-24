function [passes, updatedModel, updatedModelT, updatedBlockedMets, updatedViableMets, solvedRxns] = testSolution_type2(model, modelT, blockedMets, checkMets, viableMets, testRxns, demandRxns)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
updatedModel = model;
updatedModelT = modelT;
updatedBlockedMets = blockedMets;
updatedViableMets = viableMets;
updatedCheckMets = checkMets;
passes = false;
solvedRxns = {};

for i = 1:length(testRxns(:,1))
    mets2add = string([]);
    stoic2add = [];
    rxn2add = string(testRxns(i,1));
    compartment = "[" + extractBetween(rxn2add,"[","]") + "]";
    demandRxns_comp = demandRxns(find(contains(demandRxns,compartment)));
    updatedBlockedMets_comp = updatedBlockedMets(find(contains(updatedBlockedMets,compartment)));
    updatedViableMets_comp = updatedViableMets(find(contains(updatedViableMets,compartment)));
    updatedCheckMets_comp = updatedCheckMets(find(contains(updatedCheckMets,compartment)));
    addMet = string(testRxns{i,2}).';
    mets2add = [mets2add;addMet];
    for j = 1:length(addMet)
        stoic2add = [stoic2add;-1];
    end
    addMet = string(testRxns{i,3}).';
    mets2add = [mets2add;addMet];
    for j = 1:length(addMet)
        stoic2add = [stoic2add;1];
    end
    if cell2mat(testRxns(i,5)) == 1
        rev = true;
        lb = -1000;
        ub = 1000;
    else
        rev = false;
        lb = 0;
        ub = 1000;
    end
    
    rxn2add = convertStringsToChars(rxn2add);
    mets2add = convertStringsToChars(mets2add);
    
    [tempModelT] = addReaction(updatedModelT, rxn2add,'metaboliteList',...
                              mets2add, 'stoichCoeffList',stoic2add,...
                              'reversible', rev, 'lowerBound',lb,...
                              'upperBound', ub);

                          
    tempBlockedMets_comp = findBlockedDemandReaction(tempModelT, demandRxns_comp);
    
    disp("%%%%%%%%%%%%%%%%%%");
    disp("Initial Compartmental Blocked Count:");
    disp(length(updatedBlockedMets_comp));
    disp("Updated Compartmental Blocked Count:");
    disp(length(tempBlockedMets_comp));
    disp("%%%%%%%%%%%%%%%%%%");
    
    if length(tempBlockedMets_comp) < length(updatedBlockedMets_comp)
        disp("PASS");
        updatedModelT = tempModelT;
        updatedModel = addReaction(updatedModel, rxn2add,'metaboliteList',...
                              mets2add, 'stoichCoeffList',stoic2add,...
                              'reversible', rev, 'lowerBound',lb,...
                              'upperBound', ub);

        for k = 1:length(updatedBlockedMets_comp)
            met = updatedBlockedMets_comp(k);
            if ~any(strcmp(met,tempBlockedMets_comp))                    
                updatedViableMets = [updatedViableMets;met];
                updatedCheckMets = updatedCheckMets(updatedCheckMets~=met);
                updatedBlockedMets = updatedBlockedMets(updatedBlockedMets~=met);

            end
        end
        solvedRxns = vertcat(solvedRxns,testRxns(i,:));
        passes = true;

    else
        disp("FAIL");
    end
end
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
%       - Brandt Bessell 09/20 - Added demand reaction exclusitivity

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