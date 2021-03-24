function [blockedMets,viableMets] = findBlockedMetabolites_par(model, numWorkers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pull out intitial lists of blocked/unblocked metabolites by adding a 
%demand reaction and attempting to set as objective function for a model 
%where all media uptake reactions are allowed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelT = allowAllUptake(model);
testMets = modelT.mets;

rc =rem(length(testMets),numWorkers);
num2remove = rc;
while rc > 0 
    testMets = [testMets;"cpd99999[c]"];
    rc = rc-1;
end

checkDistSize = floor(length(testMets)/numWorkers);
checkMets_temp = strings(numWorkers,checkDistSize);

for W = 1:numWorkers
    %Divide into seperate pools for parallel processing
    if W == 1
        testMets_temp(W,:) = testMets(1:(checkDistSize*W));
    %elseif W == numWorkers
        %checkMets_temp(W,:) = checkMets((checkDistSize*(W-1)):end);
    else
        testMets_temp(W,:) = testMets((checkDistSize*(W-1)+1):(checkDistSize*W));
    end
    
blockedMets_1 = [];
blockedMets_2 = [];
blockedMets_3 = [];
blockedMets_4 = [];
viableMets_1 = [];
viableMets_2 = [];
viableMets_3 = [];
viableMets_4 = [];



parfor i = 1:length(testMets)
    via_temp = []
    blocked_temp = []
    mets = testMets
    tmet = testMets(i)
    model_temp = modelT;
    try
        model_temp = addDemandReaction(model_temp,tmet);
        %model_temp = changeObjective(model_temp, ("DM_"+ string(tmet)));
    catch 
        warning("Could not add demand reaction, would likely result in duplication");
    end
%     sol = optimizeCbModel(model_temp,'max');
%     if sol.f > 0
%         via_temp = [via_temp;tmet];
%     elseif sol.f == 0
%         blocked_temp = [blocked_temp;tmet];
%     else 
%         warning("Could not classify metabolite");
%     end

    if W == 1
        blockedMets_1 = blocked_temp;
        viableMets_1 = via_temp;
    elseif W == 2
        blockedMets_2 = blocked_temp;
        viableMets_2 = via_temp;
    elseif W == 3
        blockedMets_3 = blocked_temp;
        viableMets_3 = via_temp;
    elseif W == 4
        blockedMets_4 = blocked_temp
        viableMets_4 = via_temp;
    end
    
end
blockedMets = [blockedMets_1;blockedMets_2;blockedMets_3;blockedMets_4];
viableMets = [viableMets_1;viableMets_2;viableMets_3;viableMets_4];
blockedMets = string(sort(blockedMets));
viableMets = string(sort(viableMets));
end
