function [blockedMets,viableMets, modelT] = findBlockedMetabolites_old(model, addDemand)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pull out intitial lists of blocked/unblocked metabolites by adding a 
%demand reaction and attempting to set as objective function for a model 
%where all media uptake reactions are allowed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if addDemand
    modelT = allowAllUptake(model);
else
    modelT = model;
end
testMets = modelT.mets;
blockedMets = [];
viableMets = [];
for i = 1:length(testMets)
    tmet = testMets(i);
    disp(tmet)
    try
        if addDemand
            modelT = addDemandReaction(modelT,tmet);
        end
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
viableMets = string(sort(viableMets));
toc
end
