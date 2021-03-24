function [blockedMets,viableMets] = findBlockedMetabolites_benchmark(model, blockedMets_unsolvable, viableMets_solved)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pull out intitial lists of blocked/unblocked metabolites by adding a 
%demand reaction and attempting to set as objective function for a model 
%where all media uptake reactions are allowed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelT = allowAllUptake(model);
testMets = modelT.mets;
disp(class(testMets))
testMets = cellstr(setdiff(testMets, blockedMets_unsolvable));
disp(class(testMets))
blockedMets = [];
viableMets = [];
for i = 1:length(testMets)
    tmet = testMets(i);
    try
        modelT = addDemandReaction(modelT,tmet);
        modelT = changeObjective(modelT, ("DM_"+ string(tmet)));
    catch 
        warning("Could not add demand reaction, would likely result in duplication");
    end
    sol = optimizeCbModel(modelT,'max');
    if sol.f > 0
        viableMets = [viableMets;tmet];
    elseif sol.f == 0
        blockedMets = [blockedMets;tmet];
    else 
        warning("Could not classify metabolite");
    end
end
blockedMets = string(sort(blockedMets));
%viableMets = string(sort(viableMets));
viableMets = setdiff(viableMets_solved, blockedMets);
end
