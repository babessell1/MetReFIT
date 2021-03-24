function [passes, updatedModel, updatedmodelT, updatedBlockedMets, updatedViableMets, updatedCheckMets] = testSolution_type1(model, modelT, blockedMets, viableMets, checkMets, addID, addRxn, addRev, demandRxns)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
updatedModel = model;
updatedModelT = modelT;
updatedBlockedMets = blockedMets;
updatedViableMets = viableMets;
updatedCheckMets = checkMets;
passes = false;

tempModelT = updatedModelT;
compartment = "[" + extractBetween(addRxn(1),"[","]") + "]"; % later can make to handle multiple compartments for transports
demandRxns_comp = demandRxns(find(contains(demandRxns,compartment)));
updatedBlockedMets_comp = updatedBlockedMets(find(contains(updatedBlockedMets,compartment)));
updatedViableMets_comp = updatedViableMets(find(contains(updatedViableMets,compartment)));
updatedCheckMets_comp = updatedCheckMets(find(contains(updatedCheckMets,compartment)));

for i = 1:length(addID)
    mets2add = string([]);
    stoic2add = zeros([]);
    for j = 1:length(addRxn)
        stoic2add = [];
        spl = strsplit(addRxn(j));
        rxnside = 0;
        for k = 1:length(spl)
            try
                stoic = str2double(spl(k));
                if ~rxnSide
                    stoic2add = [stoic2add;-stoic];
                else
                    stoic2add = [stoic2add;stoic];
                end
            catch
                if spl(k) == "<=>" || spl(k) == "=>" || spl(k) == "->"
                    rxnside = 1;
                else
                    mets2add = [mets2add;spl(k)];
                end
            end
        end
    end
    
    if addRev(i) == "<=>"
        rev = true;
        lb = -1000;
        ub = 1000;
    else
        rev = false;
        lb = 0;
        ub = 1000;
    end
    [tempModelT] = addReaction(tempModelT, rxn2add,'metaboliteList',...
                              mets2add, 'stoichCoeffList',stoic2add,...
                              'reversible', rev, 'lowerBound',lb,...
                              'upperBound', ub);
end
tempBlockedMets_comp = findBlockedDemandReaction(tempModelT, 'ibm_cplex', demandRxns_comp);

if length(tempBlockedMets_comp) < length(updatedBlockedMets_comp)
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
    passes = true;
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

