function [addRxns,blockedMets,viableMets,model, modelT] = metREFit_all(model,modelT,presentDB,superDB,orthoDB,transDB,blockedMets,viableMets,demandRxns,masterMetsDB)

dbstop if error
savedWorkspace = 'MetRefit_saveData_addByRank.mat';
%reqConf = userParameters.reqConf;              %normal confidence required to pass

S_matrix = model.S;
metabolites = model.mets;
%reactions = model.rxns;
met_status = zeros(1,length(metabolites));
met_production_cnt = met_status;
%rxn_status = zeros(1,length(reactions));
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

clearvars priority_1_mets priority_2_mets priority_3_mets priority_4_mets priority_5_mets
clearvars sorted_metabolites sorted_met_prod_cnt indices met_status metabolites S_matrix 

%presRxns = presentDB.Abbreviation;
% presID = presentDB.Abbreviation;
% masterRxnIDs = string(masterDB.Abbreviation);
% masterRxns = (masterDB.Description);
% masterRxnIDsOld = string(masterDB.OldAbbreviation);
masterMetIDs = string(masterMetsDB.Abbreviation);
masterMetIDsOld = (masterMetsDB.old_abbrev);
masterMets = (masterMetsDB.NewDescription);
masterMetsOld = (masterMetsDB.Description);


% outFile = fopen('Results\TracerOutput_addByRank.txt','w');
excludeRxns = string([]);
escapeFlag = 0;
addRxns = {};
tic

%REGULAR
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
        testRxns_ortho = {};
        testRxns_rank = {};
        [testRxns_temp] = addFromSuper(met, traceNum, viableMets,...
                                       excludeRxns,superDB,orthoDB,presentDB);
                                       
        if isempty(testRxns_temp)
            fprintf("No reactions producing %s found in reference model\n", met);
            continue;
        end 
        
        for r = 1:length(testRxns_temp(:,1))
            cscore = cell2mat(testRxns_temp(r,9));
            missReacts = cell2mat(testRxns_temp(r,10));
            
            if cscore < 1 && missReacts == 0
                testRxns_ortho = vertcat(testRxns_rank,testRxns_temp(r,:));
            elseif missReacts == 0
                testRxns_rank = vertcat(testRxns_rank,testRxns_temp(r,:));              
            end
        end
        
        if isempty(testRxns_rank) && isempty(testRxns_ortho)
            fprintf("No reactions within confidence range producing %s found in reference model\n", met);
            continue;
        end
        
        solvedFlag = 0;
        
        if ~isempty(testRxns_ortho)
            [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
                = testSolution_type2(model, modelT, blockedMets, checkMets,...
                                     viableMets, testRxns_ortho, demandRxns);
        end
         
        if ~solvedFlag && ~isempty(testRxns_rank)
                [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
                    = testSolution_type3(model, modelT, blockedMets,...
                        checkMets, viableMets, testRxns_rank, demandRxns);

        end
                 
        if solvedFlag
            addRxns = vertcat(addRxns, solvedRxns);
            save(savedWorkspace, 'model', 'modelT', 'addRxns');
            numSol = length(solvedRxns(:,1));
            fprintf("%Found rank validated solution for %s\n\n", numSol, met);
            disp(solvedRxns(:,1));

        else
            fprintf("No solutions for %s could be validated\n\n", met);
        end
            
    end  
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
        fprintf("No more solutions could be found using ranked priority from reference model")
    end
end


%TRANSPORT
% while ~escapeFlag
%     initialBlockedCnt = length(blockedMets);
%     for m = 1:length(checkMets)
%         if ~any(strcmp(checkMets(m),blockedMets))
%             continue
%         end
% 
%         initMet = checkMets(m);
%         [initMetName] = retrieveMetName(initMet, masterMets,...
%                             masterMetsOld, masterMetIDs, masterMetIDsOld);
% 
%         if iscell(initMetName)
%             initMetName = initMetName{1};
%         end
%         if initMet == initMetName
%             initMetName = [];
%         end
%         
%         fprintf("\r\rAttempting to solve for metabolite, %s\n\n", initMetName);
%         met = initMet;
%         testRxns_ortho = {};
%         testRxns_rank = {};
%         [testRxns_temp] = addFromSuper(met, traceNum, viableMets,...
%                                        excludeRxns,superDB,orthoDB,presentDB);
%                                        
%         if isempty(testRxns_temp)
%             fprintf("No reactions producing %s found in reference model\n", met);
%             continue;
%         end 
%         
%         for r = 1:length(testRxns_temp(:,1))
%             cscore = cell2mat(testRxns_temp(r,9));
%             missReacts = cell2mat(testRxns_temp(r,10));
%             
%             if cscore < 1 && missReacts == 0
%                 testRxns_ortho = vertcat(testRxns_rank,testRxns_temp(r,:));
%             elseif missReacts == 0
%                 testRxns_rank = vertcat(testRxns_rank,testRxns_temp(r,:));              
%             end
%         end
%         
%         if isempty(testRxns_rank) && isempty(testRxns_ortho)
%             fprintf("No reactions within confidence range producing %s found in reference model\n", met);
%             [testRxns_temp] = addTransporter(met,traceNum,viableMets,excludeRxns,transDB,presentDB);
%             if ~isempty(testRxns_temp)
%                 beep
%                 [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
%                     = testSolution_type3(model, modelT, blockedMets,...
%                         checkMets, viableMets, testRxns_temp, demandRxns);
%                 if ~solvedFlag
%                     continue;
%                 else
%                     addRxns = vertcat(addRxns, solvedRxns);
%                     save(savedWorkspace, 'model', 'modelT', 'addRxns');
%                 end
%             else
%                 continue;
%             end
%         end
%         
%         solvedFlag = 0;
%         
%         if ~isempty(testRxns_ortho)
%             [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
%                 = testSolution_type2(model, modelT, blockedMets, checkMets,...
%                                      viableMets, testRxns_ortho, demandRxns);
%         end
%          
%         if ~solvedFlag && ~isempty(testRxns_rank)
%                 [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
%                     = testSolution_type3(model, modelT, blockedMets,...
%                         checkMets, viableMets, testRxns_rank, demandRxns);
% 
%         end
%         
%         if ~solvedFlag
%             [testRxns_temp] = addTransporter(met,traceNum,viableMets,excludeRxns,transDB,presentDB);
%             if ~isempty(testRxns_temp)
%                 [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
%                     = testSolution_type3(model, modelT, blockedMets,...
%                         checkMets, viableMets, testRxns_temp, demandRxns); 
%             end
%         end
%         
%         if solvedFlag
%             addRxns = vertcat(addRxns, solvedRxns);
%             save(savedWorkspace, 'model', 'modelT', 'addRxns');
%             numSol = length(solvedRxns(:,1));
%             fprintf("%Found rank validated solution for %s\n\n", numSol, met);
%             disp(solvedRxns(:,1));
% 
%         else
%             fprintf("No solutions for %s could be validated\n\n", met);
%         end
%             
%     end  
%     fprintf("Refreshing blocked Metabolite List\n\n")
%     [blockedMets,viableMets, modelT, demandRxns] = findBlockedMetabolites(model);
%     finalBlockedCnt = length(blockedMets);
%     if finalBlockedCnt < initialBlockedCnt
%         disp("%%%%%%%%%%%%%%%%%%")
%         disp("%%%%%%%%%%%%%%%%%%")
%         disp("Initial Blocked Count:")
%         disp(initialBlockedCnt)
%         disp("Updated Blocked Count:")
%         disp(initialBlockedCnt)
%         disp("%%%%%%%%%%%%%%%%%%")
%         disp("%%%%%%%%%%%%%%%%%%")
%         fprintf("Attempting to fill in gaps using updated blocked lists...")
%         checkMets = setdiff(checkMets_restore,viableMets);
%         escapeFlag = 0;
%     else            
%         escapeFlag = 1;
%         fprintf("No more solutions could be found using ranked priority from reference model")
%     end
% end








%%LONGCHAIN
% escapeFlag = 0;
% while ~escapeFlag
%     initialBlockedCnt = length(blockedMets);
%     for m = 1:length(checkMets)
%         if ~any(strcmp(checkMets(m),blockedMets))
%             continue
%         end
% 
%         initMet = checkMets(m);
%         [initMetName] = retrieveMetName(initMet, masterMets,...
%                             masterMetsOld, masterMetIDs, masterMetIDsOld);
% 
%         if iscell(initMetName)
%             initMetName = initMetName{1};
%         end
%         if initMet == initMetName
%             initMetName = [];
%         end
%         
%         fprintf("\r\rAttempting to solve for metabolite, %s\n\n", initMetName);
%         met = initMet;
%         testRxns = {};
%         rxnChain_1 = {};
%         [testRxns_temp] = addFromSuper(met, traceNum, viableMets,...
%                                        excludeRxns,superDB,orthoDB,presentDB);
%                                        
%         if isempty(testRxns_temp)
%             fprintf("No reactions producing %s found in reference model\n", met);
%             continue;
%         end       
% 
%         for r = 1:length(testRxns_temp(:,1))
%             rxnChain_1 = {};
%             rxnChain_1 = vertcat(rxnChain_1,testRxns_temp(r,:)); 
%             cscore = cell2mat(testRxns_temp(r,9));
%             missReacts = testRxns_temp{r,10};
%             excludeRxns = [testRxns_temp{r,1}];
%             chainSolve = 0;
%             solvedFlag = 0;
%             if missReacts == 1
%                 traceNum = traceNum + 1;
%                 chainMet_1 = testRxns_temp{r,17};
%                 rxnChain_2 = {};
%                 if length(chainMet_1) == 1
%                     [rxnChain_2] = addFromSuper(chainMet_1, traceNum, viableMets,...
%                                     excludeRxns,superDB,orthoDB,presentDB);
%                 end
% 
%                 if ~isempty(rxnChain_2)
%                     for s = 1:length(rxnChain_2(:,1))
%                         cscore = cell2mat(rxnChain_2(s,9));
%                         missReacts = cell2mat(rxnChain_2(s,10));
%                         if missReacts == 0
%                             rxnChain_temp_1 = vertcat(rxnChain_1,rxnChain_2(s,:));
%                             [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
%                                 = testSolution_type4(model, modelT, blockedMets, checkMets,...
%                                      viableMets, rxnChain_temp_1, demandRxns);
%                              if solvedFlag
%                                  break;
%                              end
%                                      
%        
%                         
%                         elseif missReacts ==1
%                             traceNum = traceNum + 1;
%                             chainMet_2 = rxnChain_2{s,17};
% % %                             chainMet_2 = testRxns_temp{s,17}; %pulls
% % info from reations for first unsolved met, not a downstream blocked
% % metabolite
%                             rxnChain_3 = {};
%                             if length(chainMet_2) == 1
%                                    [rxnChain_3] = addFromSuper(chainMet_2, traceNum, viableMets,...
%                                                 excludeRxns,superDB,orthoDB,presentDB);
%                             end
% 
%                             if ~isempty(rxnChain_3)
%                                 for t = 1:length(rxnChain_3(:,1))
%                                     cscore = cell2mat(rxnChain_3(t,9));
%                                     missReacts = cell2mat(rxnChain_3(t,10));
%                                     if missReacts == 0
%                                         rxnChain_temp_2 = vertcat(rxnChain_1,rxnChain_2,rxnChain_3(t,:));
%                                         [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
%                                             = testSolution_type4(model, modelT, blockedMets, checkMets,...
%                                                  viableMets, rxnChain_temp_2, demandRxns);
% %                                 [rxnChain_2] = addFromSuper(chainMet_2, traceNum, viableMets,...
% %                                                 excludeRxns,superDB,orthoDB,presentDB);
% %                             end
% % 
% %                             if ~isempty(rxnChain_2)
% %                                 for t = 1:length(rxnChain_2(:,1))
% %                                     cscore = cell2mat(rxnChain_2(t,9));
% %                                     missReacts = cell2mat(rxnChain_2(t,10));
% %                                     if missReacts == 0
% %                                         rxnChain_temp_2 = vertcat(rxnChain_1,rxnChain_2,rxnChain_3(t,:));
% %                                         [solvedFlag, model, modelT, blockedMets, viableMets, solvedRxns]...
% %                                             = testSolution_type4(model, modelT, blockedMets, checkMets,...
% %                                                  viableMets, rxnChain_temp_2, demandRxns);          
%                                  
%                                     end
%                                     if solvedFlag
%                                         chainSolve = 1;
%                                         break;
%                                     end
%                                 end
%                                 
%                             end                                                 
%                         end
%                         if solvedFlag 
%                              chainSolve = 1;
%                              break;
%                         end
%                     end
%                 end
%                 if chainSolve
%                     addRxns = vertcat(addRxns,rxnChain_temp_1);
%                     
%                     break;
%                 end
%             end
%         end 
%                  
%         if chainSolve
%             addRxns = vertcat(addRxns, rxnChain_temp_1);
%             save(savedWorkspace, 'model', 'modelT', 'addRxns');
%             numSol = length(rxnChain_temp_1(:,1));
%             fprintf("%Found rank validated solution for %s\n\n", numSol, met);
%             disp(rxnChain_temp_1(:,1));
% 
%         else
%             fprintf("No solutions for %s could be validated\n\n", met);
%         end
%             
%     end  
%     fprintf("Refreshing blocked Metabolite List\n\n")
%     [blockedMets,viableMets, modelT, demandRxns] = findBlockedMetabolites(model);
%     finalBlockedCnt = length(blockedMets);
%     if finalBlockedCnt < initialBlockedCnt
%         disp("%%%%%%%%%%%%%%%%%%")
%         disp("%%%%%%%%%%%%%%%%%%")
%         disp("Initial Blocked Count:")
%         disp(initialBlockedCnt)
%         disp("Updated Blocked Count:")
%         disp(initialBlockedCnt)
%         disp("%%%%%%%%%%%%%%%%%%")
%         disp("%%%%%%%%%%%%%%%%%%")
%         fprintf("Attempting to fill in gaps using updated blocked lists...")
%         checkMets = setdiff(checkMets_restore,viableMets);
%         escapeFlag = 0;
%     else            
%         escapeFlag = 1;
%         fprintf("No more solutions could be found using ranked priority from reference model")
%     end
% end


end
