function [solvedRxns,blockedMets,viableMets] = metREFit_main_bare(model, modelT, presentDB,masterDB,superDB,orthoDB,blockedMets,viableMets,demandRxns,masterMetsDB)

dbstop if error
savedWorkspace = 'MetRefit_saveData.mat';
%reqConf = userParameters.reqConf;              %normal confidence required to pass

S_matrix = model.S;
metabolites = model.mets;
reactions = model.rxns;
met_status = zeros(1,length(metabolites));
met_production_cnt = met_status;
rxn_status = zeros(1,length(reactions));
traceNum = 0;

fprintf("Setting Metabolite Statuses...\n")
num_metabolites = length(metabolites);
for i =  1:num_metabolites
    met = metabolites(i);
    if any(strcmp(met, blockedMets))
        met_status(1,i) = 1;
    else
        met_status(1,i) = 0;
    end
end
% fprintf("Setting Reaction Statuses\n")
% num_rxns = length(reactions);
% for i =  1:num_rxns
%     rxn = reactions(i);
%     if any(strcmp(rxn,blockedRxns))
%         rxn_status(1,i) = 1;
%     else
%         rxn_status(1,i) = 0;
%     end
% end
fprintf("Determining Metabolite Production Counts...\n")
disp(length(S_matrix(:,2)));
for i = 1:length(metabolites)
    for j = 1:length(S_matrix(:,i))
        if S_matrix(j,i) == 1
            met_production_cnt(i) = met_production_cnt(i) + 1;
        end
    end
end
[sorted_met_prod_cnt, indices] = sort(met_production_cnt);
sorted_metabolites = metabolites(indices);
    
fprintf("Determining Metabolite Priority...\n")
priority_1_mets = string([]);
priority_2_mets =  string([]);
priority_3_mets =  string([]);
priority_4_mets =  string([]);
priority_5_mets =  string([]);
for i = 1:num_metabolites
    prodcnt = sorted_met_prod_cnt(i);
    met = metabolites(i);
    if any(strcmp(met, blockedMets))
        if prodcnt == 0
            priority_1_mets = [priority_1_mets;metabolites(i)];
        elseif prodcnt == 1
            priority_2_mets = [priority_2_mets;metabolites(i)];
        elseif prodcnt == 2
            priority_3_mets = [priority_3_mets;metabolites(i)];
        elseif prodcnt == 3
            priority_4_mets = [priority_4_mets;metabolites(i)];
        else 
            priority_5_mets = [priority_5_mets;metabolites(i)];
        end
    end
end

checkMets = [priority_1_mets;priority_2_mets;priority_3_mets;priority_4_mets;priority_5_mets];
checkMets = flipud(checkMets);
checkMets_restore = checkMets;

presRxns = presentDB.Abbreviation;
presID = presentDB.Abbreviation;
masterRxnIDs = string(masterDB.Abbreviation);
masterRxns = (masterDB.Description);
masterRxnIDsOld = string(masterDB.OldAbbreviation);
masterMetIDs = string(masterMetsDB.Abbreviation);
masterMetIDsOld = (masterMetsDB.old_abbrev);
masterMets = (masterMetsDB.NewDescription);
masterMetsOld = (masterMetsDB.Description);


outFile = fopen('Results\TracerOutput_newDB.txt','w');
excludeRxns = string([]);
escapeFlag = 0;
addRxns = {};
tic
refresh = 0;

while ~escapeFlag
    initialBlockedCnt = length(blockedMets);
    for m = 1:length(checkMets)
        if ~any(strcmp(checkMets(m),blockedMets))
            continue
        end

        initMet = checkMets(m);
        [initMetName] = retrieveMetName(initMet, masterMets,...
                            masterMetsOld, masterMetIDs, masterMetIDsOld);

        if iscell(initMetName)
            initMetName = initMetName{1};
        end
        if initMet == initMetName
            initMetName = [];
        end
        
        fprintf("\r\rAttempting to solve for metabolite, %s\n\n", initMetName);
        met = initMet;
        testRxns = {};
        [testRxns_temp] = addFromSuper(met, traceNum, viableMets,...
                                       excludeRxns,superDB,orthoDB,presentDB);
                                       
        if isempty(testRxns_temp)
            fprintf("No reactions producing %s found in reference model\n", met);
            continue;
        end 
        
        for r = 1:length(testRxns_temp(:,1))
            cscore = cell2mat(testRxns_temp(r,9));
            missReacts = cell2mat(testRxns_temp(r,10));
            
            if cscore == 1 && missReacts == 0
                testRxns = vertcat(testRxns,testRxns_temp(r,:));
            end
        end
        
        if isempty(testRxns)
            fprintf("No orthologically verified reactions producing %s found in reference model\n", met);
            continue;
        end
            
        [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
            = testSolution_type2(model, modelT, blockedMets,checkMets,...
                                 viableMets, testRxns, demandRxns);
          
        if solvedFlag
            addRxns = vertcat(addRxns, solvedRxns);
            save(savedWorkspace, 'model', 'modelT', 'addRxns');
            numSol = length(solvedRxns(:,1));
            fprintf("%d orthologically validated solutions found solving for %s\n\n", numSol, met);
            disp(solvedRxns(:,1));
            refresh = 1;
        else
            fprintf("No solutions for %s could be orthologically validated\n\n", met);
        end
            
    end
    
    if refresh
        fprintf("Refreshing blocked Metabolite List\n\n")
        [blockedMets,viableMets, modelT, demandRxns] = findBlockedMetabolites(model);
        finalBlockedCnt = length(blockedMets);
        if finalBlockedCnt < initialBlockedCnt
            disp("%%%%%%%%%%%%%%%%%%")
            disp("%%%%%%%%%%%%%%%%%%")
            disp("Initial Blocked Count:")
            disp(initialBlockedCnt)
            disp("Updated Blocked Count:")
            disp(initialBlockedCnt)
            disp("%%%%%%%%%%%%%%%%%%")
            disp("%%%%%%%%%%%%%%%%%%")
            fprintf("Attempting to fill in gaps using updated blocked lists...")
            checkMets = setdiff(checkMets_restore,viableMets);
            escapeFlag = 0;
        else            
            escapeFlag = 1;
            fprintf("No more solutions could be found using reference model")
        end
        refresh = 0;
    end
end
end
