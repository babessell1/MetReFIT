function [tracerResults,trialInfo,dataStore] = tracer2(presentDB,masterDB,superDB,orthoDB,blockedMets,viableMets,masterMetsDB,originData, userParameters, userMod)
%function [ID2add,rxn2add,rev2add,boostID2add,boostRxn2add,boostRev2add,trialInfo,dataStore, tStart, tEnd] = tracer2(presentDB,masterDB,superDB,orthoDB,blockedMets,viableMets,masterMetsDB,userParameters, userMod)
%Needed for function to work, ignore tot variables here
%totgene = 0;totX = 0;totT = 0;tot0 = 0;tot1 = 0;tot2 = 0;tot3 = 0;tot4 = 0;tot5 = 0;tot6 = 0;rcnt=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  USER-DEFINED PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%
%MAXIMUM NUMBER OF TRIALS permitted for solving each blockage (default=30)
%Lower to save time, increase for extra thoroughness (no effect on score)
%% 
dbstop if error

%%% TODO
% change fndBlockedMetabolites to accept a model where demand rxns have already been added and parallize
% If orthologic match to Zea maize --> Cscore = 1.00
% Add new Section:
%%%If solution passes with score of 1:0 (genetic match):
%%%%% make temp model (copy model)
%%%%% add solution to temp model
%%%%% check if target metabolite will now solve
%%%%%%%% if yes run findBlockedMetabs and ensure blocked list is smaller
%%%%%%%%%%%%%% if yes add solution to model and model w/ demand reactions and update blocked list
%%%%%%%%%%%%%% if no -- run again again allowing for a lower threshold score
%%%%%%%%%%%%%%%%%%5 if fails, give up and continue








stdTrials = userParameters.stdTrials;
extraTrials = userParameters.extraTrials;           %Extra Trials for rulebreaking
masterTrials = userParameters.masterTrials;         %Trials if masterdatabase is being used (may need to be very large for complete coverage of possibilities)
largePathTrials = userParameters.largePathTrials;   %Trials for large path generation 

% stdTrials = 50;
% extraTrials = 100;          %Extra Trials for rulebreaking
% masterTrials = 500;         %Trials if masterdatabase is being used (may need to be very large for complete coverage of possibilities)
% largePathTrials = 1000;     %Trials for large path generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAX PERMISSIBLE CONFIDENCE SCORE for solution to pass (default 3.0)
%Increasing will speed up solution and allow a wider range of potential 
%pathways to pass at the expense of confidence.

reqConf = userParameters.reqConf;              %normal confidence required to pass
boostConf = userParameters.boostConf;
masterConf = userParameters.mastConf;
largePathConf = userParameters.largePathConf;

% reqConf = 3.0;              %normal confidence required to pass
% boostConf = 10.0;
% mastConf = 30.0;
% largePathConf = 15.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Skipping checking of various database FUTURE FEATURE

% skAdThres = 100;
% skSuThres = 200;
% skRand = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RULEBREAKING -lets unsolved paths to try again under different, more
%lenient restrictions and/or with increased number of permutations. 
%Pentalty is added to confidence scores of solutions breaking specific rule

allowExTrials = userParameters.allowExTrials;
allowExConf = userParameters.allowExConf;
allowMaster = userParameters.allowMaster;


% easementPenalty = @(cscore) 0.5*cscore;
% extraPenalty = @(cscore) 0.1*cscore;


%LARGE-SCALE NETWORKING AND PATH LENGTH (Unstable)

%Turns on network trial length permanantly
largePathPerm = userParameters.largePathPerm;
% largePathPerm = 0;
%Turns on larger networks for unsolvable metabolites
largePathTemp = userParameters.largePathTemp;
% largePathTemp = 0;
%maximum currently unresolved metabolites per set of reactions downstream
%from target 
netRxnsMax = userParameters.netRxnsMax;
% netRxnsMax = 2;
%maximum length of pathway or network, measured as number metabolites
%downstream from initial metabolite (default = 3) *THIS IS IN EFFECT WHEN LARGE SCALE NET IS OFF
stdLenMax = userParameters.stdLenMax;
% stdLenMax = 2;
%maximim length of pathway for conditionally used network finding
largeLenMax = userParameters.largeLenMax;
% largeLenMax = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  CONFIDENCE SCORE CALCULATION  %%%%%%%%%%%%%%%%%%%%%%

%Variable Names (not exhaustive)

%cscore            confidence score (variable being calculated)

%Equation


%%%%%%%%%%%%%%%%%%%  END OF USER DEFINED PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%First Flag to append rulebroken rxns to end of list, this should be safe
%without the second flag as the components of the reactions should now all
%be synthesized in the model.

%second allows easement reactions to be added to the boosted listed.   
%*This is NOT stable, these reactions could very well contain metabolites
%that are not being synthesized by the model even with all new reactions 
%pulled from this function. They are better served for prediction purposes
%and thus is not recomended for adding without a manual vetting process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTER LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Starting from the final metabolite in the list of blocked
%metabolites run initial tracerLoop function to pull out precurser for
%first downstream reaction. If output is acceptable and all precursers are
%not already viable, admit reactants into the inner loop.

% checkMets = blockedMets; %metabolites which will be checked
% checkMets = "cpd11255_mr[p]";
if userMod == "standard"
    checkMets = blockedMets;
elseif userMod == "exampleMalate" || userMod == "testMalate"
    checkMets = "cpd00130_mr[m]"; %Testing malate synthesis from fumarate in TCA cycle
    viableMets = [viableMets;"cpd00001_mr[m]";"cpd00106_mr[m]"];
elseif userMod == "exampleAnthocyanin" || userMod == "testAnthocyanin"
    checkMets = "cpd03501_mr[c]";     %Testing anthocyanin production pathway
%     checkMets = "cpd08919_mr[c]"; %Testing anthocyanin production pathway
elseif userMod == "examplePyruvate" || userMod == "testPyruvate"  %Currently not incorporated
    checkMets = "cpd00020_mr[c]";   %Testing glycolysis pathway (removal of pyruvate and phosphoenolpyruvate production)
elseif userMod == "exampleUnsatFattyAcids" || userMod == "testUnsatFattyAcids"  %Currently not incorporated
    checkMets = ["cpd03848_mr[c]"; "cpd03848_mr[p]"; "cpd03848_mr[x]"; "cpd00188_mr[c]"; "cpd00188_mr[c]"; "cpd00327_mr[x]"; "cpd14883_mr[x]"; "cpd00070_mr[x]"; "cpd11812_mr[c]"]; %Example showing filling of biosynthesis of Unsaturated fatty acids pathway   %Testing glycolysis pathway (removal of pyruvate and phosphoenolpyruvate production)
end

presRxns = presentDB.Abbreviation;
presID = presentDB.Abbreviation;
masterRxnIDs = string(masterDB.Abbreviation);
masterRxns = (masterDB.Description);
masterRxnIDsOld = string(masterDB.OldAbbreviation);
masterMetIDs = string(masterMetsDB.Abbreviation);
masterMetIDsOld = (masterMetsDB.old_abbrev);
masterMets = (masterMetsDB.NewDescription);
masterMetsOld = (masterMetsDB.Description);
superRxns = superDB.Reactions;



if largePathPerm
    largePathFlag = 1;
    maxTrials = largePathTrials;
    largePathFlag = 1;
    maxLength = largeLenMax;
else
    largePathFlag = 0;
    maxLength = stdLenMax;
    maxTrials = stdTrials;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exTrialFlag = 0; %allows additional number of trials to be attempted
exConfFlag = 0; %allows higher confidences scores pass
masterFlag = 0; %allows master database to be utilized
outeridx = 0; %index used for choosing which reaction to us. When multiple
%reactions are found in the tracerLoop function, a large list may be found
% so the list is iterated through until a match is found that satisfies a
% confidence score, if none are found, this value will return to 0 and the
%search may be restarted with new parameters if the user settings allow for it

innerRandFlag = 0; %flag that allows for random combinations to be tested, it is only really used when the large path discovery setting is on because there are so many possible combinations -> we should look into an alternative to using a randomized system
idx_shift = 0;     %Index that indicates how many reactions to shift the inner trial list. 
avgPoss = 1;
trial = 1;  %counter for how many trials have been tested

%file handling and some output to let you know its started
%outFile = fopen('D:\\cobratoolbox\TracerOutput2020_1.txt','w');
outFile = fopen('Results\TracerOutput_newDB.txt','w');

fprintf(outFile,"TRACER-GAPFILL RESULTS\n\n");
fprintf("\r\rProcessing Blocked Metabolites...\n\n");
possList=[];
tic
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
    %!FIX
    if trial > stdTrials/2 || trial > avgPoss %rem(trial,2)==0 %skipAdd tells the tracerLoop to skip checking for existing reactions first, I think I was probably trying to force skip for even values to diversify the potential solutions but this needs be handled better. Regardless, the ability to skip checking existing reactions at the begining of the search is critical for covering all posibilities
        skipAdd=1;
    else
        skipAdd=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract initial reactions, iterate through until end until suitible 
    %solution is found if exhausted all iteretions before then, randomize
    %index for inner indicies
    [addRxns] = tracerLoop2(met,traceNum,tempViaMets,tempExcludeRxns,presentDB,superDB,orthoDB,masterDB,largePathFlag,masterFlag,netRxnsMax,skipAdd,presRxns,blockedMets);
%     warning("check here")
    if ~iscell(addRxns)
        if addRxns == "SKIP"
           if ~isempty(initMetName)
               fprintf("RESULTS FOR METABOLITE %s\t%s:\n\rAll requirements are met for reaction, reason flux is not simulating is unknown\n\n",initMet,initMetName);
           else
               fprintf("RESULTS FOR METABOLITE %s\t:\n\rAll requirements are met for reaction, reason flux is not simulating is unknown\n\n",initMet);
           end
            escapeFlag = 1;
            skipFlag=1;
        elseif ~largePathFlag && largePathTemp
            warning("FATALERROR");
            largePathFlag = 1;
            outeridx = 0;
            maxTrials=largePathTrials;
        elseif largePathFlag 
            warning("FATALERROR");
            escapeFlag = 1;
        else
            fprintf("FAILUREUNKNOWN");
            escapeFlag=1;
        end
    elseif isempty(addRxns)
        if ~isempty(initMetName)
            fprintf("RESULTS FOR METABOLITE %s\t%s:\n\r",initMet,initMetName);
        else
            fprintf("RESULTS FOR METABOLITE %s\t:\n\r",initMet);
        end
        
        if~largePathFlag && largePathTemp
            fprintf("No downstream reactions found under standard requirements, reattempting with less stringent parameters\n");
            largePathFlag = 1;
            maxLength=largeLenMax;
            outeridx = 0;
        elseif largePathFlag
            fprintf("No downstream reactions found under easement, temporarily disabling metabolite from the blocked list\nAdding metabolite to seperate output list, consider manual validation of associated reaction(s)\n");
            escapeFlag = 1;
        else
            fprintf("No downstream reactions found under standard requirements,passing\n");
            escapeFlag=1;
        end
    else
        nextLen=length(addRxns(:,1)); %number of qualifying reactions found
        possib=nextLen; %number of possibilities to test
        if nextLen>outeridx %iterates through all possible solutions 
            outeridx=outeridx+1;
        else
            %!FIX
            outeridx = 1; %after trying all potential solutions for the first reaction, the selected first solution is randomized -> I do not like this approach an alternative should be found here
            innerRandFlag = 1;
            idx_shift = idx_shift-1;
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
                outerComps = outerComps + "(" + prodt(i) + "*)";
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
            tempBlockMets = tempBlockMets(~contains(tempBlockMets,prodt));
            tempViaMets = unique([tempViaMets;prodt]);
            rulelist(end)= 'x';
            solvedFlag = 1;
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
        while ~isempty(loopMets) && ~forceStop && pathLength<=maxLength
            lastLen=nextLen; 
            traceNum = traceNum+1;
            newLoopMets = string([]);
            for i = 1:length(loopMets)
                met = loopMets(i);
                if met ~= "X"
                    rulelist(end) = 'n';
                end
                [addRxns] = tracerLoop2(met,traceNum,tempViaMets,tempExcludeRxns,presentDB,superDB,orthoDB,masterDB,largePathFlag,masterFlag,netRxnsMax,skipAdd,presRxns);
                
                %move reactions who have products that are both missing
                %metabolites to the top of the list
%                 if ~isempty(addRxns) && length(addRxns(:,1))>1 && ~skipAdd
%                     for i = 1:length(addRxns(:,1))
%                         rt = addRxns{i,2}; %reactants
%                         pt = addRxns{i,3}; %products
%                         rx = addRxns{i,1}; %reaction id
%                         if all(ismember(pt,blockedMets)
                
                %Randomize Index after nonrandomized attempts fail            
                if iscell(addRxns) && ~isempty(addRxns) && pathLength<maxLength
                    pcheck = char(addRxns{1,6}); %checks the origin of the first/best qualifying reaction from tracerLoop
                    inTrialList=[];
                    for j=1:length(addRxns(:,1))
                        inTrialList=[inTrialList;j];
                    end

                    if innerRandFlag %%%&& contains(pcheck,'0')%%% %if the randomizer flag from before is active (meaning no solution was found for all of the outerloop designated starting point reaction) AND the first/best qualifying reations for the next reaction is a reaction already in the model,
                        %idx_shift = idx_shift - 1;
                        inTrialList = circshift(inTrialList, idx_shift);
                    end
                    
                    chidx=1;
                    passFlag=0;
                    nextLen=length(inTrialList);
                    possib= nextLen*lastLen; %total number of possibilities
                    while chidx<=nextLen && ~passFlag
                            innertrial=inTrialList(chidx);
                            chidx=chidx+1;
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
                    forceStop = 1;
                    rulelist(end+1) = 'f';
                end   
            end
            checkMets = flipud(checkMets);
            checkMets = flipud([checkMets;"X"]);
            loopMets = newLoopMets;
            if isempty(loopMets) && ~forceStop && rlist(end)~=rlist(1)
                rulelist(end+1) = 'x';
                %Can be turned on to play a noise when a solution was found
                %I originally used it for testing, can be removed in final
                %build
%                 if rlist(end)==rlist(1)
%                      load handel
%                      sound(y,Fs)
%                    warning('stop');
%                 end
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
%                 NR = 9;  % max rank, will need to change to a user setting later
%                 Np = 5; %will need to change this to be a user setting later
%                 R = 0;    %reaction rank from superDB
%                 sym = string(extractBetween(tlist(i),"(",")"));
%                 if contains(sym, "+") %GPR match
%                     R = 9;
%                     %B = 1 so,
%                     nG = str2double(erase(sym,"+"));
%                     gFactor = nG/Np;
%                     pFactor = R/NR;   %P
%                 elseif sym == "T" || sym == "X" %Transfer reaction or Master reaction
%                     %B = 0 so,
%                     gFactor = 0;
%                     pFactor = 2;   %P
%                 elseif sym == "0" %Present model rxn
%                     gFactor = 0;
%                     pFactor = 0;
%                     modTotalRxns = modTotalRxns + 1;
%                 elseif sym == "1" || sym == "2" || sym == "3" || sym == "4" || sym == "5"|| sym == "6" || sym == "7" || sym == "8" || sym == "9"
%                     R = str2num(sym);
%                     gFactor = nG/Np;
%                     pFactor = R/NR;   %P
%                 else
%                     error("INVALID TRACER SYMBOL\nSymbols must be formated such that the traceNum is directly between parenthesis (ex. '<(#)>' ) \n\n");
%                 end
%                 C(i) = (pFactor/4 + 3*gFactor/4);
            end
            cscore = sum(cscore_total)/(totalRxns-modTotalRxns); %totalRxns = T
            
            

           
            
%             if exTrialFlag && ~easementFlag
%                 cscore = cscore + extraPenalty(cscore);
%             elseif easementFlag
%                 cscore = cscore + easementPenalty(cscore);
%             end
    
            %confidence scores that do not pass have the option of being
            %resolved under different parameters specified in the user
            %config at the begining of the script
            try
                %if ( ((cscore >= reqConf) && ~exConfFlag) || ((cscore >=boostConf) && exConfFlag && ~masterFlag) || ((cscore>=masterConf) && masterFlag) || ((cscore>=largePathConf) && largePathFlag) ) && (~all(contains(rlist,presRxns))) || skipAdd
                if cscore >= reqConf
                    if ~exTrialFlag && ~exConfFlag && ~masterFlag && ~largePathFlag 
                        fprintf(outFile,"\rRESULTS FOR METABOLITE %s\t%s\n\nSolved under standard parameters\n\n", initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved under standard parameters\n\n",initMet,initMetName);
                    elseif exTrialFlag && ~exConfFlag
                        fprintf(outFile, "\rRESULTS FOR METABOLITE %s\t%s\n\nSolved given additional permutation\n\n", initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved given additional permutation\n\n", initMet,initMetName);
                    elseif exConfFlag && ~masterFlag
                        fprintf(outFile, "\rRESULTS FOR METABOLITE %s\t%s\n\nSolved for increased confidence threshold\n\n",initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved for increased confidence threshold\n\n",initMet,initMetName);
                    elseif masterFlag && ~largePathFlag
                        fprintf(outFile,"\rRESULTS FOR METABOLITE %s\t%s\n\nSolved using Master model\n\n",initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved using Master model\n\n",initMet,initMetName);
                    elseif largePathFlag

                    else
                        %error("\rRESULTS FOR METABOLITE "+initMet+"\n\rSolved under unreasonable parameters\n\r");
                    end
                    escapeFlag = 1;
                    solvedFlag = 1;
                end
            catch
                warning("Unknown Problem Causing this if statement to fail")
            end
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
                    exTrialFlag = 1;
                end             
            elseif skipAdd && trial>2*avgPoss
                fprintf(trial +"\n");
%                  avgPoss=mean(possList);
                fprintf("\naverage # of Possibilities: %s\n",string(avgPoss))
                trial=999999;
                exTrialFlag = 1;
            end
            if trial>=maxTrials 
                fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nNo Solutions could be determined in the alotted number of trials,\n",initMet,initMetName);
                if ~allowExTrials && ~allowExConf && ~allowMaster && ~largePathTemp
                    fprintf("Standard (%d) was insufficient for a complete solution. Consider enabling rulebreaking and/or manual validation of associated reaction(s)\n\n",stdTrials);
                    escapeFlag=1;
                    idx_shift = 0;
                elseif allowExTrials && ~exTrialFlag
                    fprintf("Standard (%d) was insufficient for a complete solution.\nReattempting with extra trials\n\rTrial:\n",stdTrials);
                    exTrialFlag = 1;
                    maxTrials = extraTrials;
                    trial = 1;
                    outeridx=0;
                    idx_shift = 0;
                elseif allowExConf && ~exConfFlag
                    if trial ~= 999999
                        fprintf("Increasing max number of trials was still insufficient for a complete solution.\nReattempting with less permissible confidence\n\rTrial:\n");
                    else
                        fprintf("Standard (%d) was insufficient for a complete solution.\nExtra trials will also fail.\nReattempting with less permissible confidence\n\rTrial:\n",stdTrials);
                    end
                    exConfFlag = 1;
                    trial = 1;
                    outeridx=0;    
                    idx_shift = 0;
                elseif allowMaster && ~masterFlag
                    fprintf("Increasing permissible confidence was still insufficient for a complete solution.\nReattempting with use of master database\n\rTrial:\n");
                    masterFlag = 1;
                    trial = 1;
                    outeridx=0;
                    idx_shift = 0;
                elseif largePathTemp && ~largePathFlag
                    fprintf("Utilization of master list was still insufficient for a complete solution.\nReattempting for larger and/or more open pathways\n\rTrial:\n");
                    outeridx=0;
                    trial=1;
                    largePathFlag=1;
                    maxTrials = largePathTrials;
                    idx_shift = 0;
                elseif largePathTemp && largePathFlag
                    fprintf("Extra trials and easement of parameters was still insufficient for a complete solution\nTemporarily disabling metabolite from the blocked list\nAdding metabolite to seperate output list,\nconsider manual validation of associated reaction(s)\n\n ");
                    escapeFlag = 1;               
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
       innerRandFlag = 0;
       trialStore = trial;
       trial = 1;
       idx_shift = 0; 
       if solvedFlag && ~all(contains(rlist,presRxns))
           viableMets = tempViaMets;
           trialSuccess = 1; 
                      
           if ~masterFlag && ~largePathFlag
               blockedMets = tempBlockMets;
               checkMets(~contains(checkMets,tempBlockMets)) = [];
               if exConfFlag
                   boostMets = [boostMets;initMet];
                   boostID2add = [boostID2add;tempAddID];
                   boostRxn2add = [boostRxn2add;tempAddRxn];
                   boostRev2add = [boostRev2add;tempAddRev];
                   Btr2add = [Btr2add;tempTrace];
               else
                   passedMets = [passedMets;initMet];
                   rxn2add = [rxn2add;tempAddRxn];
                   ID2add = [ID2add;tempAddID];
                   rev2add = [rev2add;tempAddRev];
                   trace2add = [trace2add;tempTrace];
               end
           elseif masterFlag && ~largePathFlag
               mastMets = [mastMets;initmet];
               checkMets = checkMets(1:end-1);
           elseif largePathFlag
               largePathMets = [largePathMets;initMet];
           end
           %Placeholder Output
           
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
                   %printTrace = "(" + originData(str2double(extractBetween(tlist(j),'(',')'))).origin + ")";
%                elseif tlist(j) == "(1)"
%                    printTrace = "(Reaction originated from the model species (super database))";
%                elseif tlist(j) == "(2)"
%                    printTrace = "(Reaction originated from species of rank 2 (super database))";
%                elseif tlist(j) == "(3)"
%                    printTrace = "(Reaction originated from species of rank 3 (super database))";
%                elseif tlist(j) == "(4)"
%                    printTrace = "(Reaction originated from species of rank 4 (super database))";
%                elseif tlist(j) == "(5)"
%                    printTrace = "(Reaction originated from species of rank 5 (super database))";
%                elseif tlist(j) == "(6)"
%                    printTrace = "(Reaction originated from species of rank 6 (super database))";
               end
                    
%                    printTrace = tlist(j);
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
       exTrialFlag=0;
       exConfFlag=0;
       masterFlag=0;
       if ~largePathPerm
            maxTrials = stdTrials;
            largePathFlag=0;
            pathLength = stdLenMax;
       else
           maxTrials = largePathTrials;
        end
    end
    checkMets = checkMets(~strcmp(checkMets,"X"));
    
    fullRxnIDs = [outerRxnIDs+innerRxnIDs];
    fullTracers = [outerTracers+innerTracers];
    fullRxns = [outerRxns+innerRxns];
    fullRevs = [outerRevs+innerRevs];
    fullComps = [outerComps+innerComps];
    fullTraceNums = [outerTraceNums+innerTraceNums];
    
    %%%%WARNING: CHANGE PLACEMENT TO BEFORE FLAGS ARE RESET--[RESOLVED]
    
%     if innerRandFlag
%         randPermActive = 1;
%     end
%     if skipAdd
%         skipAddActive = 1;
%     end
    
    

end
%Write to Output file
goodScore = [transpose(ID2add);transpose(trace2add);transpose(rxn2add);transpose(string(rev2add))];
fprintf(outFile,"\n\n\nReactions solved with Confidence Score under 3.0 (Better)\r\n");
fprintf(outFile,'%s %s %s %s\r\n', 'ID','Tracer','Rev','Reaction');
fprintf(outFile,'"%s" "%s" "%s" "s\r\n',goodScore);

badScore = [transpose(boostID2add);transpose(Btr2add);transpose(boostRev2add);transpose(boostRxn2add)];
fprintf(outFile,"\n\n\nReactions solved with Confidence Score between 3.0 and 5.0\r\n");
fprintf(outFile,'%s %s %s %s\r\n', 'ID','Tracer','Rev','Reaction');
fprintf(outFile,'"%s" "%s" "%s" "s\r\n',badScore);

fprintf(outFile,"\r\nPassed Mets\r\n");
fprintf(outFile, '%s', passedMets);
fprintf(outFile,"\r\nBoosted Mets\r\n");
fprintf(outFile, '%s', boostMets);
fprintf(outFile,"\r\nFailed Mets\r\n");
fprintf(outFile, '%s', failedMets);
fprintf(outFile,"\r\nSkipped Mets\r\n");
fprintf(outFile, '%s', skippedMets);
fclose(outFile);
warning('STOP');

tracerResults.ID2add = ID2add;
tracerResults.rxn2add = rxn2add;
tracerResults.rev2add = rev2add;
tracerResults.boostID2add = boostID2add;
tracerResults.boostRxn2add = boostRxn2add;
tracerResults.boostRev2add = boostRev2add;
toc
end

        



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    