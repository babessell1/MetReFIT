function [tracerResults,trialInfo,dataStore] = metREFit_main(model, modelT, presentDB,masterDB,superDB,orthoDB,blockedMets,viableMets,demandRxns, masterMetsDB,originData, userParameters)

dbstop if error
type1cnt = 0;
type2cnt = 0;


allowMaster = userParameters.allowMaster;
allowTransport = 0;
stdTrials = userParameters.stdTrials;
extraTrials = userParameters.extraTrials;           %Extra Trials for rulebreaking
masterTrials = userParameters.masterTrials;         %Trials if masterdatabase is being used (may need to be very large for complete coverage of possibilities)

reqConf = userParameters.reqConf;              %normal confidence required to pass

S_matrix = model.S;
metabolites = model.mets;
reactions = model.rxns;
met_status = zeros(1,length(metabolites));
met_production_cnt = met_status;
rxn_status = zeros(1,length(reactions));

fprintf("Setting Metabolite Statuses\n")
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
fprintf("Determining Metabolite Production Counts\n")
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

%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE OUTPUT LISTS  %%%%%%%%%%%%%%%%%%%%%%%%%%
ID2add = [];
rxn2add = [];
rev2add = [];
boostID2add = [];
boostRxn2add = [];
boostRev2add = [];
trace2add = [];
Btr2add = [];
boostMets =[];
failedMets = [];
passedMets = [];
skippedMets = [];

presRxns = presentDB.Abbreviation;
presID = presentDB.Abbreviation;
masterRxnIDs = string(masterDB.Abbreviation);
masterRxns = (masterDB.Description);
masterRxnIDsOld = string(masterDB.OldAbbreviation);
masterMetIDs = string(masterMetsDB.Abbreviation);
masterMetIDsOld = (masterMetsDB.old_abbrev);
masterMets = (masterMetsDB.NewDescription);
masterMetsOld = (masterMetsDB.Description);

%Starting from the final metabolite in the list of blocked
%metabolites run initial tracerLoop function to pull out precurser for
%first downstream reaction. If output is acceptable and all precursers are
%not already viable, admit reactants into the inner loop.
stdLenMax = userParameters.stdLenMax;
maxLength = stdLenMax;
maxTrials = stdTrials;
masterFlag = 0; %allows master database to be utilized
outeridx = 0; %index used for choosing which reaction to us. When multiple

idx_shift = 0;     %Index that indicates how many reactions to shift the inner trial list. 
avgPoss = 1;
trial = 1;  %counter for how many trials have been tested

outFile = fopen('Results\TracerOutput_newDB.txt','w');
fprintf("\r\rProcessing Blocked Metabolites...\n\n");
possList=[];
tic
%[outRxns] = tracerExistReactExtract(met,presentDB,traceNum);

trySuper_outer = 0;
tryTransport_outer = 0;
trySuper_inner = 0;

while ~isempty(checkMets)
    %Assigning value and initializing values
    cscore  = 0;        
    pathLength = 0;                 %stores the path length which is equal to the number of downstream reactions so far added to find a solution
    escapeFlag = 0;                 %used to force exit the outer loop
    solvedFlag = 0;                 %tells the loop to exit and also save and display the found solution
    skipFlag = 0;                   %flags when the reason for flux simulation not working is unknown because all required metabolites are presumably viable, exits the loop
    forceStop = 0;                  %this might not do anything anymore 
    traceNum = 1;                   %number used to track where reactions and metabolites lie downstream from the initial metabolite
    initMet = checkMets(end);       %initial metabolite whose blockage is being solved for
    initMetName = [];
    [initMetName] = retrieveMetName(initMet, masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
    if iscell(initMetName)
        initMetName = initMetName{1};
    end
    if initMet == initMetName
        initMetName = [];
    end
    
    met = initMet; %met to be passed through additional functions, starts as the initial but will change as the solver goes downstream
    tempBlockMets = blockedMets;
    %these temporary lists provide a way for the true lists to update as
    %new reactions and metabolites are solved for and added. Once a set of
    %reactions is passed as a solution these lists are used to update the
    %master list, otherwise they are discarded
    tempViaMets = viableMets;
    tempExcludeRxns = string([]);
    tempAddID = [];
    tempAddRev = [];
    tempAddRxn = [];
    tempTrace = [];
    currentMets = [];
    %these list are all used for formating the visualzed output pathway
    tlist = []; %contains tracer symbol values for reactions being added
    clist = []; %contains compound values
    rlist = []; %contains reaction ID values
    rklist =[]; %contains trace number values
    revlist = []; %contains reversiblity values
    rxlist = []; %contains reaction values
    rankList = [];
    rulelist = ''; %list of characters that tell the function how to format: c = present compound, s = target/blocked compound, n = new reaction, x = terminator
    conflist = [];
    testRxns = {};
    testFlag = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract initial reactions, iterate through until end until suitible 
    %solution is found if exhausted all iterations before then, randomize
    %index for inner indicies

    if ~trySuper_outer && ~tryTransport_outer 
        [addRxns] = addFromExisting(met, traceNum, tempViaMets, tempExcludeRxns, presentDB);
    elseif ~tryTransport_outer
        [addRxns] = addFromSuper(met,traceNum,tempViaMets,tempExcludeRxns,superDB,orthoDB,presentDB);
        if ~testFlag && ~isempty(addRxns)
            for i = 1:length(addRxns(:,1))
                if addRxns{i,10} == 0 && addRxns{i,9} == 1.0
                    testRxns = vertcat(testRxns,addRxns(i,:));              
                end
                testFlag = 1;
            end
        elseif testFlag  && ~isempty(addRxns)
            for i = 1:length(addRxns(:,1))
                if addRxns{i,10} == 0 && addRxns{i,9} == 1.0
                addRxns{i,:} = [];
                end
                testFlag = 0;
            end
        end
    else
        [addRxns] = addTransport(met, traceNum, tempViaMets, tempExcludeRxns);
    end

      
    if ~iscell(addRxns)
        if addRxns == "SKIP"
           if ~isempty(initMetName)
               fprintf("RESULTS FOR METABOLITE %s\t%s:\n\rAll requirements are met for reaction, reason flux is not simulating is unknown\n\n",initMet,initMetName);
           else
               fprintf("RESULTS FOR METABOLITE %s\t:\n\rAll requirements are met for reaction, reason flux is not simulating is unknown\n\n",initMet);
           end
           escapeFlag = 1;
           skipFlag=1;
        else
            fprintf("FAILED");
            escapeFlag=1;
        end
        
    elseif isempty(addRxns)
        if ~trySuper_outer
            trySuper_outer = 1;
        elseif ~tryTransport_outer && allowTransport
            tryTransport_outer = 1 ;
        else
            escapeFlag=1;
            fprintf("No solutions found under standard requirements,passing\n");
        end
        
    else
        outerLen=length(addRxns(:,1)); %number of qualifying reactions found
        possib=outerLen; %number of possibilities to test
        if outerLen>outeridx %iterates through all possible solutions 
            outeridx=outeridx+1;
        elseif ~trySuper_outer
            trySuper_outer = 1;
            forceStop = 1;
        elseif ~tryTransport_outer
            forceStop = 1;
        end
        
        %after a reaction is selected, its information is extracted into
        %the followingl lists which are used in the inner loop where the
        %next reaction will be selected using a similar algorithm
        rxnID = addRxns{outeridx,1}; 
        react = transpose(addRxns{outeridx,2});
        prodt = transpose(addRxns{outeridx,3});
        rxn = addRxns{outeridx,4};
        rev = addRxns{outeridx,5};
        trace = addRxns{outeridx,6};
        confscore = addRxns{outeridx,9};
        trNumber = addRxns{outeridx,7}; 
        %counters are reset to 0
        ccnt=0;
        scnt=0;
        rcnt=0;
        loopMets = []; %metabolites which pass solution requirements in the inner loop will populate this list, an empty list shows that no solution was found and the loop should run again using updated parameters
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Initialize rulelist, a character based code for translating the
        %pathways into a visualized format
        for i=1:length(prodt)
     
            currentMets=[currentMets;prodt(i)];        
            if prodt(i) == initMet
                [metName] = retrieveMetName(prodt(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                clist = [clist;"*"+metName+"*"];
                rulelist(end+1) = 'c';
                ccnt=ccnt+1;
                tempViaMets=[tempViaMets;prodt(i)];
            elseif ismember(prodt(i),tempViaMets)
                [metName] = retrieveMetName(prodt(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                clist = [clist;metName];
                rulelist(end+1) = 's';
                ccnt = ccnt+1;
            else
                [metName] = retrieveMetName(prodt(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                clist = [clist;metName+"*"];
                %outerComps = outerComps + "(" + prodt(i) + "*)";
                rulelist(end+1) = 'c';
                scnt=scnt+1;
                tempViaMets=[tempViaMets;prodt(i)];
            end
        end
        rulelist(end+1) = 'n';
        for i=1:length(react)
            currentMets=[currentMets;react(i)];
            if ismember(react(i),viableMets)
                [metName] = retrieveMetName(react(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                clist = [clist;metName];
                metIdx = find(contains(masterMetIDs, react(i)));
                rulelist(end+1) = 's';
                ccnt = ccnt+1;
            else
                [metName] = retrieveMetName(react(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                clist = [clist;metName+"*"];
                rulelist(end+1) = 'c';
                loopMets = [loopMets;react(i)];
                scnt=scnt+1;
                tempViaMets=[tempViaMets;react(i)];
            end
        end
        rulelist(end+1) = 'n';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Update lists of reaction information for visual translation,
        %store information not already contained in the model in temporary
        %lists to be used if a solution is found
        tlist = [tlist;trace];
        conflist = [conflist;confscore];
        rklist = [rklist;trNumber];
        revlist = [revlist;rev];  
        rxlist = [rxlist;rxn];
        tempExcludeRxns = [tempExcludeRxns;rxnID];   
                
        if ~any(strcmp(presID,rxnID))
            tempAddID = [tempAddID;rxnID];
            tempAddRxn = [tempAddRxn;rxn];
            tempAddRev = [tempAddRev;rev];
            tempTrace = [tempTrace;trace];
            rlist = [rlist;rxnID+"*"];
            rcnt=rcnt+1;
        else
            rlist = [rlist;rxnID];
        end
           
        if isempty(loopMets)
            if tlist(1) ~= "(0)"
                tempBlockMets = tempBlockMets(~contains(tempBlockMets,prodt));
                tempViaMets = unique([tempViaMets;prodt]);
                rulelist(end)= 'x';
                solvedFlag = 1;
            else 
                trySuper_outer = 1;
            end
        else
            forceStop = 0;
            tempBlockMets = tempBlockMets(~contains(tempBlockMets,prodt));
            tempViaMets = unique([tempViaMets;prodt]);
        end
        pathLength = 1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        %%%%%%%%%%%%%%%%%%%%%%%  INNER LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       repeat process of outer loop for reactants of initial metabolite,
%       repeat until the tracerloop function fails to find a suitable
%       reaction or until the reactants pulled from the function are
%       already viable. *Note that the loop operates off of 3 seperate
%       lists of blocked metabolites. The blockedMetabolites is preserves
%       the integrity of the list and is ONLY updated if the final
%       reactants in the path downstream of the inital metabolite are 
%       viable. The temp list allows for the list to be updated and used
%       within the loop so that metabolites are not being inaccurately 
%       referenced under the assumption that the downstream reactions are
%       being added. The checkMets list is the least exclusive, being 
%       updated even in the event that a metabolite is not able to be 
%       solved, this allows the function to progress without getting hung
%       up on an unsolvable metabolite. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        while ~isempty(loopMets) && ~forceStop && pathLength<=maxLength && ~testFlag
            traceNum = traceNum+1;
            newLoopMets = string([]);
            for i = 1:length(loopMets)
                met = loopMets(i);
                if met ~= "X"
                    rulelist(end) = 'n';
                end
                  
                if ~trySuper_outer && ~tryTransport_outer 
                    [addRxns] = addFromExisting(met, traceNum, tempViaMets, tempExcludeRxns, presentDB);
                elseif ~tryTransport_outer
                    [addRxns] = addFromSuper(met,traceNum,tempViaMets,tempExcludeRxns,superDB,orthoDB,presentDB); 
                elseif allowTransport
                    [addRxns] = addTransport(met, traceNum, tempViaMets, tempExcludeRxns);
                end
                
                if iscell(addRxns) && ~isempty(addRxns) && pathLength<maxLength
                    pcheck = char(addRxns{1,6}); %checks the origin of the first/best qualifying reaction from tracerLoop
                    inTrialList=[];
                    for j=1:length(addRxns(:,1))
                        inTrialList=[inTrialList;j];
                    end
                    
                    inneridx=1;
                    passFlag=0;
                    innerLen=length(inTrialList);
                    possib= innerLen*outerLen; %total number of possibilities
                    while inneridx<=innerLen && ~passFlag
                            innertrial=inTrialList(inneridx);
                            inneridx=inneridx+1;
                            rxn = addRxns{innertrial,4};
                            react = transpose(addRxns{innertrial,2});
                            prodt = transpose(addRxns{innertrial,3});
                            rxnID = addRxns{innertrial,1};
                        if (~all(ismember(prodt,currentMets)) && ~all(ismember(react,currentMets))) || contains(rxnID, "Transport")%loop passes with reaction information if products and reactants are not all currently viable, this ensures that the function does not get stuck using a reversible reaction as the solution to the previous set of needed metabolites forcing the into an infinite loop
                            rxnID = addRxns{innertrial,1};
                            react = transpose(addRxns{innertrial,2});
                            prodt = transpose(addRxns{innertrial,3});
                            rxn = addRxns{innertrial,4};
                            rev = addRxns{innertrial,5};
                            trace = addRxns{innertrial,6};
                            trNumber = addRxns{innertrial,7};
                            passFlag=1;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Update lists for visual translation
                    if passFlag
                        for j=1:length(react)
                            currentMets=[currentMets;react(j)];                      
                            if  ismember(react(j),viableMets)
                                [metName] = retrieveMetName(react(j), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                                clist = [clist;metName];
                                rulelist(end+1) = 's';
                                scnt = scnt+1;
                            else
                                [metName] = retrieveMetName(react(j), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                                clist = [clist;metName+"*"];
                                metIdx = find(contains(masterMetIDs, react(j)));
                                newLoopMets = [newLoopMets;react(j)];
                                rulelist(end+1) = 'c';
                                ccnt = ccnt+1;
                                tempViaMets=[tempViaMets;react(j)];
                            end
                        end
                        clist = [clist;"X"];
                        tlist = [tlist;trace];
                        conflist = [conflist;confscore];
                        revlist = [revlist;rev];
                        rxlist = [rxlist;rxn];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        tempExcludeRxns = [tempExcludeRxns;rxnID];
                        if ~any(strcmp(presID,rxnID))
                            tempAddID = [tempAddID;rxnID];
                            tempAddRxn = [tempAddRxn;rxn];
                            tempAddRev = [tempAddRev;rev];
                            tempTrace = [tempTrace;trace];
                            rlist = [rlist;rxnID+"*"];
                            rcnt = rcnt+1;
                        else
                            rlist = [rlist;rxnID];
                        end
                    end
                else
                    if ~trySuper_inner
                        trySuper_inner = 1;
                        trySuper_outer = 0;
                        tryTransport_outer = 0;
                    else
                        allowTransport = 0;
                    end
                    forceStop = 1;
                    rulelist(end+1) = 'f';
                end   
            end
            checkMets = flipud(checkMets);
            checkMets = flipud([checkMets;"X"]);
            loopMets = newLoopMets;
            if isempty(loopMets) && ~forceStop && rlist(end)~=rlist(1)
                rulelist(end+1) = 'x';

            end
            pathLength = pathLength+1;

        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %If a pathway is sucessfully solved to completion, calculate its 
        %confidence score. If acceptable solution is found, escape the loop,
        %reset starting parameters,commit changes to blocked/viable lists 
        %and add rxns + corresponding IDs to output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Summation of tracer values
        if rulelist(end) == 'x'%         totgene = 0;
            cscore_total = 0;
            totalRxns = length(tlist);
            modTotalRxns = 0;
            for i = 1:totalRxns % calc individual scores for each added reaction
                cscore_total = cscore_total + conflist(i);
                sym = string(extractBetween(tlist(i),"(",")"));
                if sym == "0"
                     modTotalRxns = modTotalRxns + 1;
                end
            end
            cscore = sum(cscore_total)/(totalRxns-modTotalRxns); %totalRxns = T
            %confidence scores that do not pass have the option of being
            %resolved under different parameters specified in the user
            %config at the begining of the script
            if cscore >= reqConf && testFlag
                [passesFlag, model, modelT, blockedMets, viableMets, checkMets] = testSolution_type2(model, modelT, blockedMets, checkMets, viableMets, testRxns, demandRxns);
                type = 2;
                testRxns = {};
                
            elseif cscore >= reqConf
                type = 1;
                [passesFlag, model, modelT, blockedMets, viableMets, checkMets] = testSolution_type1(model, modelT, blockedMets, viableMets, checkMets, addID, addRxn, addRev, demandRxns);
            else
                testRxns = {};
                
            end
                

            try
                if cscore >= reqConf && passesFlag
                    
                    if ~masterFlag
                        fprintf(outFile,"\rRESULTS FOR METABOLITE %s\t%s\n\nSolved under standard parameters\n\n", initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved under standard parameters\n\n",initMet,initMetName);
                    else
                        fprintf(outFile,"\rRESULTS FOR METABOLITE %s\t%s\n\nSolved using Master model\n\n",initMet,initMetName);
                    end
                    escapeFlag = 1;
                    solvedFlag = 1;
                end
            catch
                warning("Unknown Problem Causing this if statement to fail")
            end
        elseif testFlag 
            
        end
        %Advance to next trial if failed to produce an acceptable
        %pathway. If max trial is hit,either change parameters or pass
        %metabolite
        if ~escapeFlag
            possList=[possList,possib];
            trial = trial+1;
            avgPoss=mean(possList);
            if mod(trial,10) == 0
                fprintf(trial +"\n");
%                  avgPoss=mean(possList);
                fprintf("\naverage # of Possibilities: %s\n",string(avgPoss))
                %this uses the total number of calculated possibilities
                %(averaged across all trials) to determine when to force
                %quit. 2 times the calculated possibilities are used since
                %randomization is used so this should be updated when that
                %is fixed. 
                if trial>2*avgPoss
                    trial=999999;
                end             
            elseif trial>2*avgPoss
                fprintf(trial +"\n");
%                  avgPoss=mean(possList);
                fprintf("\naverage # of Possibilities: %s\n",string(avgPoss))
                trial=999999;
            end
            if trial>=maxTrials 
                fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nNo Solutions could be determined in the alotted number of trials,\n",initMet,initMetName);
                if allowMaster && ~masterFlag
                    fprintf("Increasing permissible confidence was still insufficient for a complete solution.\nReattempting with use of master database\n\rTrial:\n");
                    masterFlag = 1;
                    trial = 1;
                    outeridx=0;
                    idx_shift = 0;             
                else
                    fprintf("Standard (%d) was insufficient for a complete solution. Consider enabling rulebreaking and/or manual validation of associated reaction(s)\n\n",stdTrials);
                    escapeFlag = 1;
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Update output information and initial loop parameters
    if escapeFlag
       possList=[];
       outeridx = 0;
       trialStore = trial;
       trial = 1;
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if solvedFlag && type == 1 && ~all(contains(rlist,presRxns))
           beep
           trialSuccess = 1; 
           type1cnt = type1cnt+1;           
           if ~masterFlag
               passedMets = [passedMets;initMet];
               rxn2add = [rxn2add;tempAddRxn];
               ID2add = [ID2add;tempAddID];
               rev2add = [rev2add;tempAddRev];
               trace2add = [trace2add;tempTrace];
           else
               mastMets = [mastMets;initmet];
           end

           gpr_idx = 0;
           for j = 1:length(rxlist)
               printID = rlist(j);
               currentRxn = rlist(j);
               strippedRxn = strip(currentRxn,"*");
               matchIdx = find(contains(masterRxnIDs,strippedRxn));
               %If rxn ID is not found in new IDs, try with old IDs
               if isempty(matchIdx)
                   matchIdx = find(contains(masterRxnIDsOld, strippedRxn));
               end
               %If ID is found, assign description.
               if ~isempty(matchIdx)
                   if ~ismissing(masterRxns(matchIdx))
                       rxnName = masterRxns(matchIdx);
                   else
                       rxnName = currentRxn;
                   end
               else
                   rxnName = currentRxn;
               end
               rlist(j) = rxnName;
               printRxn = rxlist(j);
               if contains(tlist(j),"+")
                   printTrace = "(Has GPR)";
               elseif tlist(j) == "(X)"
                   printTrace = "(Reaction originated from master database)";
               elseif tlist(j) == "(T)"
                   printTrace = "(Reaction is a transfer reaction)";
               elseif tlist(j) == "(0)"
                   printTrace = "(Reaction originated from current model)";
               elseif str2double(extractBetween(tlist(j),'(',')')) <= length(originData)
                   printTrace =tlist(j);
               end
                    
               fprintf(printID+"   "+printRxn+"   "+printTrace+"\n");
               fprintf(outFile,printID+"   "+printRxn+"   "+printTrace+"\r\n"); 
               if tlist(j) == "(+)"
                   gpr_idx = gpr_idx + 1;
                   fprintf("\tReaction is catalyzed by the following genes:\n")
                   fprintf(outFile, "\tReaction is catalyzed by the following genes:\r\n");
                   orthologGenesStored = gpr_match_origins{gpr_idx,1};                        %Retrieving appropriate genes in reference species
                   outputOrthologGene = extractBetween(orthologGenesStored,'[',']');                   %Extracting from storage format
                   targetSpeciesOrthologsStored = gpr_match_origins{gpr_idx, 2};                       %Retrieving corresponding genes in target species
                   targetSpeciesOrthologs = extractBetween(targetSpeciesOrthologsStored, '[',']');                %Extracting from storage format
                   gene_origins = gpr_match_origins{gpr_idx, 3};                    %Retrieving the origin species names from genes in reference species
                   gene_origins_extracted = extractBetween(gene_origins, '[',']');  %Extracting from storage format
                   %Displaying result
                   for gene_idx = 1:length(outputOrthologGene)
                       fprintf('\t -- %s corresponding to the model species'' gene(s) %s\n', outputOrthologGene{gene_idx}, targetSpeciesOrthologs{gene_idx})
                       fprintf(outFile,'\t -- %s corresponding to the model species'' gene(s) %s\n', outputOrthologGene{gene_idx}, targetSpeciesOrthologs{gene_idx});

                       
%                        fprintf('\t -- %s in %s corresponding to the %s gene(s) %s\n', outputOrthologGene{gene_idx}, gene_origins_extracted{gene_idx}, target_species, targetSpeciesOrthologs{gene_idx})
%                        fprintf(outFile,'\t -- %s in %s corresponding to the %s gene(s) %s\n', outputOrthologGene{gene_idx}, gene_origins_extracted{gene_idx}, target_species, targetSpeciesOrthologs{gene_idx});
                   end
               end
           end
           strscore = string(cscore);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %format output
           strcode = string(rulelist);
           fprintf("\nConfidence Score: %s\n",strscore);
           fprintf(outFile,"\r\nConfidence Score: %s\r\n",strscore);
           fprintf("\nConversion Code: %s\n",strcode);
           fprintf(outFile,"\r\nConversion Code: %s\r\n",strcode);
           fprintf(outFile,"........................................................................................................\r\n\r\n");
           fprintf("........................................................................................................\r\n");
           lineList = strsplit(strcode,"n");
           blockL = 0;
           blockR = 100;
           clistcnt=0;
           rlistcnt=0;
           llen=length(lineList);
           try
               for i = 1:length(lineList)
                   groupList = strsplit(lineList(i),"t");
                   pos=[];
                   compLen=[];
                   posCode='';
                   finishFlag=0;
                   if length(groupList)>0
                       %groupSpacing = (blockR-blockL)/length(groupList)
                       groupSpacing = (blockR-blockL)/(length(groupList));
                       groupspace=blockL;
                       for j=1:length(groupList)
                           group = groupList(j);
                           gchr=char(groupList(j));
                           if endsWith(gchr,'x')
                               finishFlag=1;
                               gchr=gchr(1:end-1);
                           end
                           groupMetLength=strlength(gchr);
                           %groupspace=groupspace+groupSpacing;
                           groupspace=(j-1)*groupSpacing;
                           for k=1:groupMetLength
                               metpos=groupspace+(groupSpacing/(groupMetLength+1));
                               pos=[pos,k*metpos/groupMetLength];
                               posCode=[posCode,gchr(k)];
                           end
                           posCode=[posCode,'t'];
                           pos=[pos,0];
                       end
                   else
                       group=gchr(i);
                       groupspace=(blockR-blockL)/2;
                       for k=1:length(group)
                           pos=[pos,k*metpos];
                           posCode=[posCode,char(group(k))];
                       end
                   end
                   for j=1:strlength(posCode)
                       if j==1
                           clistcnt=clistcnt+1;
                           compid=clist(clistcnt);
                           len=strlength(compid);
                           compLen=[compLen,len];
                           asize=round(pos(j));
                           printForm=strjoin(repmat(" ",1,asize),"")+clist(clistcnt);
                       elseif posCode(j-1)=='t'
                           clistcnt=clistcnt+1;
                           compid=clist(clistcnt);
                           len=strlength(compid);
                           compLen=[compLen,len];
                           asize=round(pos(j)-pos(j-1));
                           printForm=strjoin(repmat(" ",1,asize),"")+clist(clistcnt);
                       elseif posCode(j)~='t'
                           clistcnt=clistcnt+1;
                           compid=clist(clistcnt);
                           len=strlength(compid);
                           compLen=[compLen,len];
                           asize=round(pos(j)-pos(j-1));
                           printForm=strjoin(repmat(".",1,asize),"")+clist(clistcnt);
                       else
                           compLen=[compLen,0];
                           printForm="";
                       end
                       fprintf(printForm);
                       fprintf(outFile,printForm);
                   end
                   fprintf("\n\n")
                   fprintf(outFile,"\r\n\r\n");
                   for m=1:5
                       nextRxn=1;
                       if m==1
                           for j=1:strlength(posCode)
                               if posCode(j)=='c' && nextRxn
                                   rlistcnt=rlistcnt+1;
                               elseif posCode(j)=='t'
                                   nextRxn=1;
                               end
                               if j==1
                                   asize=round(pos(j)+(compLen(j)));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                        if revlist(rlistcnt)==1
                                            printForm=strjoin(repmat(" ",1,asize),"")+"/\\"+strjoin(repmat(" ",1,bsize),"");
                                        else
                                            printForm=strjoin(repmat(" ",1,asize),"")+"/\\"+strjoin(repmat(" ",1,bsize),"");
                                        end
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"")+"  "+strjoin(repmat(" ",1,bsize),"");
                                   end
                               else
                                   asize=round(pos(j)-pos(j-1)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                       if revlist(rlistcnt)==1
                                            printForm=strjoin(repmat(" ",1,asize),"")+"/\\"+strjoin(repmat(" ",1,bsize),"");
                                       else
                                           printForm=strjoin(repmat(" ",1,asize),"")+"/\\"+strjoin(repmat(" ",1,bsize),"");
                                       end
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"")+"  "+strjoin(repmat(" ",1,bsize),"");
                                   end
                               end
                               fprintf(printForm);
                               fprintf(outFile,printForm);
                           end
                           fprintf("\n")
                           fprintf(outFile,"\r\n");
                       elseif m==3
                           for j=1:strlength(posCode)
                               if posCode(j)=='t'
                                   nextRxn=1;
                               end
                               if j==1
                                   asize=round(pos(j)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                        printForm=strjoin(repmat(" ",1,asize),"")+rlist(rlistcnt);
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"");
                                   end
                               else
                                   asize=round(pos(j)-pos(j-1)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                        printForm=strjoin(repmat(" ",1,asize),"")+rlist(rlistcnt);
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"");
                                   end
                               end
                               fprintf(printForm);
                               fprintf(outFile,printForm);
                           end
                           fprintf("\n")
                           fprintf(outFile,"\r\n");
                       elseif m==5
                           for j=1:strlength(posCode)
                               if posCode(j)=='t'
                                   nextRxn=1;
                               end
                               if j==1
                                   asize=round(pos(j)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                        if revlist(rlistcnt)==1
                                            printForm=strjoin(repmat(" ",1,asize),"")+"\\/"+strjoin(repmat(" ",1,bsize),"");
                                        else
                                            printForm=strjoin(repmat(" ",1,asize),"")+"| "+strjoin(repmat(" ",1,bsize),"");
                                        end
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"")+"  "+strjoin(repmat(" ",1,bsize),"");
                                   end
                               else
                                   asize=round(pos(j)-pos(j-1)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                       if revlist(rlistcnt)==1
                                            printForm=strjoin(repmat(" ",1,asize),"")+"\\/"+strjoin(repmat(" ",1,bsize),"");
                                       else
                                           printForm=strjoin(repmat(" ",1,asize),"")+"| "+strjoin(repmat(" ",1,bsize),"");
                                       end
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"")+" "+strjoin(repmat(" ",1,bsize),"");
                                   end
                               end
                               fprintf(printForm);
                               fprintf(outFile,printForm);
                           end                   
                       else
                           for j=1:strlength(posCode)
                               if j==1
                                   asize=round(pos(j)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                        if revlist(rlistcnt)==1
                                            printForm=strjoin(repmat(" ",1,asize),"")+"||"+strjoin(repmat(" ",1,bsize),"");
                                        else
                                            printForm=strjoin(repmat(" ",1,asize),"")+"| "+strjoin(repmat(" ",1,bsize),"");
                                        end
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"")+"  "+strjoin(repmat(" ",1,bsize),"");
                                   end
                               else
                                   asize=round(pos(j)-pos(j-1)+compLen(j));
                                   asize=asize-round(0.5*compLen(j));
                                   bsize=round(0.5*(compLen(j)));
                                   if posCode(j)=='c' && nextRxn
                                       nextRxn=0;
                                       if revlist(rlistcnt)==1
                                            printForm=strjoin(repmat(" ",1,asize),"")+"||"+strjoin(repmat(" ",1,bsize),"");
                                       else
                                           printForm=strjoin(repmat(" ",1,asize),"")+"| "+strjoin(repmat(" ",1,bsize),"");
                                       end
                                   else
                                       printForm=strjoin(repmat(" ",1,asize),"")+" "+strjoin(repmat(" ",1,bsize),"");
                                   end
                               end
                               fprintf(printForm);
                               fprintf(outFile,printForm);
                           end
                           fprintf("\n")
                           fprintf(outFile,"\r\n");
                       end     
                   end 
                   fprintf("\n\n")
                   fprintf(outFile,"\r\n\r\n");
               end
               fprintf(outFile,"........................................................................................................\r\n\r\n");
               fprintf("........................................................................................................\r\n");    
               fprintf(outFile,"########################################################################################################\r\n\r\n");
               fprintf("########################################################################################################\r\n");
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           catch
               warning("Problem Converting into visual format")
           end
       elseif solvedFlag && type == 2
           fprintf("Solved using type 2")
           beep
           beep
           type2cnt = type2cnt+1;
           
       elseif skipFlag
           trialSuccess = 2; 
           checkMets = checkMets(1:end-1);
           skippedMets = [skippedMets;initMet];
       else
           trialSuccess = 3; 
           checkMets = checkMets(1:end-1);
           fprintf("No solution found for " + initMet + " skipping over and continuing to next\n\n");
           failedMets = [failedMets;initMet];
           
            
            
       end
       masterFlag=0;
    end
    checkMets = checkMets(~strcmp(checkMets,"X"));
    
%     fullRxnIDs = [outerRxnIDs+innerRxnIDs];
%     fullTracers = [outerTracers+innerTracers];
%     fullRxns = [outerRxns+innerRxns];
%     fullRevs = [outerRevs+innerRevs];
%     fullComps = [outerComps+innerComps];
%     fullTraceNums = [outerTraceNums+innerTraceNums];
    

    
    

end
end
