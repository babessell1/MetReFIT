function [tracerResults,trialInfo,dataStore] = tracer2_par(presentDB,masterDB,superDB,orthoDB,blockedMets,viableMets,masterMetsDB,originData, userParameters, userMod, numWorkers)
%function [ID2add,rxn2add,rev2add,boostID2add,boostRxn2add,boostRev2add,trialInfo,dataStore, tStart, tEnd] = tracer2(presentDB,masterDB,superDB,orthoDB,blockedMets,viableMets,masterMetsDB,userParameters, userMod)
%Needed for function to work, ignore tot variables here
%totgene = 0;totX = 0;totT = 0;tot0 = 0;tot1 = 0;tot2 = 0;tot3 = 0;tot4 = 0;tot5 = 0;tot6 = 0;rcnt=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  USER-DEFINED PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%
%MAXIMUM NUMBER OF TRIALS permitted for solving each blockage (default=30)
%Lower to save time, increase for extra thoroughness (no effect on score)
%% 
%Assigning model species
target_species = "zma";

%Creating timer for function
tStart = tic;

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
mastConf = userParameters.mastConf;
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

% allowExTrials = 1;
% allowExConf = 1;
% allowMaster = 0;
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
%length(rlist)     total number of reactions in the pathway
%rcnt              number of NEW reactions in the pathway
%scnt              number of already viable compounds in path
%ccnt              number of blocked compounds in path
%totgene           number of rxns with organism GPR (Zma) 
%totX              number of rxns pulled from master model 
%totT              number of transport rxns added
%tot0              number of reactions already present in model
%tot1              number of rank 1 reactions from super model(Zea Mays) 
%tot2              number of rank 2 reactions (SorgumBicolor)
%tot3              number of rank 3 reactions (Foxtail Millet)
%tot4              number of rank 4 reactions (Arabidopsis)
%tot5              number of rank 5 reactions (Brachypodium)
%tot6              number of rank 6 reactions (Aegilops, barley, Oryza) 
            
%TRACER SCORES (score given to the source of an added reaction)

% %Variable name      Definition/Source      Tracer Tags
% scoreGene = 0;      %Has GPR                  "(+)"
% scoreX = 10.0;      %Master                   "(X)"
% scoreT = 1.0;       %Transport                "(T)"
% score0 = 0;         %Current model            "(0)"
% score1 = 0;         %zma (super)              "(1)"
% score2 = 0.5;       %sbi                      "(2)"
% score3 = 0.75;      %sita                     "(3)"
% score4 = 1.0;       %aly,ath                  "(4)"
% score5 = 1.25;      %bdi                      "(5)"
% score6 = 2.0;       %obr,osa,ogl,ats,bar      "(6)"

%Equation
%default: scoreCalc = @() SUM(scoreGene*totgene+scoreX*totX+scoreT*totT
    %+score0*tot0+score1*tot1+score2*tot2+score3*tot3+score4*tot4
    %+score5*tot5+score6*tot6)*(1+0.1*rcnt);

% scoreCalc = @(totgene,totX,totT,tot0,tot1,tot2,tot3,tot4,tot5,tot6,rcnt,ccnt)...
%     sum(scoreGene*totgene+scoreX*totX+scoreT*totT...
%     +score0*tot0+score1*tot1+score2*tot2+score3*tot3+score4*tot4...
%     +score5*tot5+score6*tot6);

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



%Storage array for data used to compare versions of tracer2
dataStore = {}; 
%Creating headers for output file

%Initalizing tracker variables

randCnt = 0;                    %Tracks number of times randPerm is called w/out skipAdd active
skipRandCnt = 0;                %Tracks number of times randPerm is called w/ skipAdd active
skipCnt = 0;                    %Tracks number of times skipAdd is called w/out randPerm active
randSkipCnt = 0;                %Tracks number of times skipAdd is called w/ randPerm active
preRandRxnCntTot = 0;           %Tracks number of rxns attempted w/out randPerm or skipAdd for met
randRxnCntTot = 0;              %Tracks number of rxns attempted w/ randPerm active
skipRxnCntTot = 0;              %Tracks number of rxns attempted w/ skipAdd active
randSkipRxnCntTot = 0;          %Tracks number of rxns attempted w/ randPerm and skipAdd active
innerPreRandRxnCntTot = 0;      %Tracks number of rxns attempted from inner loop w/out randPerm and skipAdd
innerRandRxnCntTot = 0;         %Tracks number of rxns attempted from inner loop w/ randPerm
innerSkipRxnCntTot = 0;         %Tracks number of rxns attempted from inner loop w/ skipAdd
innerRandSkipRxnCntTot = 0;     %Tracks number of rxns attempted from inner loop w/ randPerm and skipAdd
outerCScoreTot = 0;             %Tracks total contribution of outer loop to confidence score of all attempted solutions
randOuterCScoreTot = 0;         %Tracks total contribution of outer loop to confidence score of all attempted solutions w/ randPerm
skipOuterCScoreTot = 0;         %Tracks total contribution of outer loop to confidence score of all attempted solutions w/ skipAdd
randSkipOuterCScoreTot = 0;     %Tracks total contribution of outer loop to confidence score of all attempted solutions w/ randPerm and skipAdd
preRandCScoreTot = 0;           %Total confidence score of all attempted solutions without randPerm and skipAdd
randCScoreTot = 0;              %Total confidence score of all attempted solutions w/ randPerm
skipCScoreTot = 0;              %Total confidence score of all attempted solutions w/ skipAdd
randSkipCScoreTot = 0;          %Total confidence score of all attempted solutions w/ randPerm and skipAdd
innerLoopCntTot = 0;            %Tracks total number of times through inner loop for met
innerForLoopCntTot  = 0;        %Total number of times through FOR loop in inner loop per met
outerLoopCnt = 0;               %Total number of times through outer loop for met


%Stores information per trial
trialInfo = {};               

for W = 1:numWorkers
    trialInfo(W, 1, 1:25) = {'Initial Metabolite' 'IDs of rxns added' 'IDs of rxns added from outer loop' 'IDs of rxns added from inner loop'...
        'Tracer symbols of rxns' 'Tracer Symbols for outer loop rxns' 'Tracer symbols for inner loop rxns'...
        'Rxns values' 'Outer rxn values' 'Inner rxn values' 'All rxns looped through in inner loop'...
        'Reversibilities' 'Reversibilities from outer loop' 'Reversibilities from inner loop'...
        'Compounds' 'Compounds from outer loop' 'Compounds from inner loop'...
        'Trace numbers' 'Trace numbers from outer loop' 'Trace numbers from inner loop'...
        'Level of rulebreaking applied' 'skipAdd active?' 'randPerm active?' 'Successful/unsuccessful trial' 'Failed due to conf score?'};

    dataStore(W,1,1:10) = {'Metabolite' 'Highest level of rulebreaking applied' 'Trials' 'Path Length' 'randPerm calls' 'randPerm calls w/ skipAdd active' 'skipAdd calls' 'skipAdd calls w/ randPerm active' 'rxns added w/out rand' 'rxns added w/ randPerm'};
    dataStore(W,1,11:20) = {'rxns added w/ skipAdd' 'rxns added w/ randPerm & skipAdd' 'total rxns attempted w/out randPerm & skipAdd' 'total rxns attempted w/ randPerm' 'total rxns attempted w/ skipAdd' 'total rxns attempted w/ both randPerm & skipAdd' 'rxns added in inner loop w/out rand' 'rxns added in inner loop w/ randPerm' 'rxns added in inner loop w/ skipAdd' 'rxns added in inner loop w/ randPerm & skipAdd'};
    dataStore(W,1,21:30) = {'total rxns attempted in inner loop w/ out rand' 'total rxns attempted in inner loop w/ randPerm' 'total rxns attempted in inner loop w/ skipAdd' 'total rxns attemped in inner loop w/ randPerm & skipAdd' 'Previously blocked mets added to solution w/out rand' 'Previously blocked mets added to solution w/ randPerm' 'Previously blocked mets added to solution w/ skipAdd' 'Previously blocked mets added to solution w/ both randPerm & skipAdd' 'Previously blocked mets (products) added to solution w/ out skipAdd' 'Previously blocked mets (products) added to solution w/ skipAdd'};
    dataStore(W,1,31:40) = {'Previously blocked mets added to solution from inner loop w/out randPerm and skipAdd' 'Previously blocked mets added to solution from inner loop w/ randPerm' 'Previously blocked mets added to solution from inner loop w/ skipAdd' 'Previously blocked mets added to solution from inner loop w/ randPerm and skipAdd' 'Confidence score of solution' 'Contribution of outer loop to confidence score of solution' 'Tot conf score contribution of outer loop' 'Tot conf score contribution of outer loop (randPerm trials)' 'Tot conf score contribution of outer loop (skipAdd trials)' 'Tot conf score contribution of outer loop (randPerm & skipAdd trials)'};
    dataStore(W,1,41:50) = {'Total conf score w/out randPerm & skipAdd' 'Total conf score w/ randPerm' 'Total conf score w/ skipAdd' 'Total conf score w/ randPerm & skipAdd' 'times through inner loop for solution' 'total times through inner loop' 'total times through FOR loop of inner loop' 'times through outer loop for metabolite' 'randPerm/skipAdd active for final trial' 'successful/unsuccesful trial'};

end    
%Initalizing tracker variables
    
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
trialFile = fopen('Results\TracerTrialInfo_newDB.txt','w');
dataFile = fopen('Results\TracerDataInfo_newDB.txt','w');
fprintf(trialFile, "Initial Metabolite\tIDs of rxns added\tIDs of rxns added from outer loop\tIDs of rxns added from inner loop\tTracer symbols of rxns\tTracer Symbols for outer loop rxns\tTracer symbols for inner loop rxns values\tFull rxn values\tOuter rxn values\tInner rxn values\tAll rxns looped through in inner loop\tReversibilities\tReversibilities from outer loop\tReversibilities from inner loop\tCompounds\tCompounds from outer loop\tCompounds from inner loop\tTrace numbers\tTrace numbers from outer loop\tTrace numbers from inner loop\tLevel of rulebreaking applied\tskipAdd active?\trandPerm active?\tSuccessful/unsuccessful trial\tConfidence score result in failure?\n");
fprintf(dataFile, 'Metabolite\tHighest Level of Rule Breaking Applied\t Trials\tPath Length\trandPerm calls\trandPerm calls w/ skipAdd active\tskipAdd calls\tskipAdd calls w/ randPerm active\trxns added w/out rand\trxns added w/ randPerm\trxns added w/ skipAdd\trxns added w/ randPerm & skipAdd\ttotal rxns attempted w/out randPerm & skipAdd\ttotal rxns attempted w/ randPerm\ttotal rxns attempted w/ skipAdd\ttotal rxns attempted w/ both randPerm & skipAdd\trxns added in inner loop w/out rand\trxns added in inner loop w/ randPerm\trxns added in inner loop w/ skipAdd\trxns added in inner loop w/ randPerm & skipAdd\ttotal rxns attempted in inner loop w/ out rand\ttotal rxns attempted in inner loop w/ randPerm\ttotal rxns attempted in inner loop w/ skipAdd\ttotal rxns attemped in inner loop w/ randPerm & skipAdd\tPreviously blocked mets added to solution w/out rand\tPreviously blocked mets added to solution w/ randPerm\tPreviously blocked mets added to solution w/ skipAdd\tPreviously blocked mets added to solution w/ both randPerm & skipAdd\tPreviously blocked mets (products) added to solution w/ out skipAdd\tPreviously blocked mets (products) added to solution w/ skipAdd\tPreviously blocked mets added to solution from inner loop w/out randPerm and skipAdd\tPreviously blocked mets added to solution from inner loop w/ randPerm\tPreviously blocked mets added to solution from inner loop w/ skipAdd\tPreviously blocked mets added to solution from inner loop w/ randPerm and skipAdd\tConfidence score of solution\tContribution of outer loop to confidence score of solution\tTot conf score contribution of outer loop\tTot conf score contribution of outer loop (randPerm trials)\tTot conf score contribution of outer loop (skipAdd trials)\tTot conf score contribution of outer loop (randPerm & skipAdd trials)\tTotal conf score w/out randPerm & skipAdd\tTotal conf score w/ randPerm\tTotal conf score w/ skipAdd\tTotal conf score w/ randPerm & skipAdd\ttimes through inner loop for solution\ttotal times through inner loop\ttotal times through FOR loop of inner loop\ttimes through outer loop for metabolite\trandPerm/skipAdd active for final trial\tsuccessful/unsuccesful trial\n'); 


fprintf(outFile,"TRACER-GAPFILL RESULTS\n\n");
fprintf("\r\rProcessing Blocked Metabolites...\n\n");

%parpool(numWorkers);
blockedMets_temp = strings(numWorkers,length(blockedMets));
viableMets_temp = strings(numWorkers,length(viableMets));
outerLoopCnt_temp = zeros(numWorkers,1);
exTrialFlag_temp = zeros(numWorkers,1);
exConfFlag_temp = zeros(numWorkers,1);
masterFlag_temp = zeros(numWorkers,1);
largePathFlag_temp = zeros(numWorkers,1);
avgPoss_temp = zeros(numWorkers,1);
trial_temp = zeros(numWorkers,1);
skipCnt_temp = zeros(numWorkers,1);
randSkipCnt_temp = zeros(numWorkers,1);
outeridx_temp = zeros(numWorkers,1);
maxTrials_temp = zeros(numWorkers,1);
maxLength_temp = zeros(numWorkers,1);
innerRandFlag_temp = zeros(numWorkers,1);
idx_shift_temp = zeros(numWorkers,1);
preRandRxnCntTot_temp = zeros(numWorkers,1);
skipRxnCntTot_temp = zeros(numWorkers,1);
outerCScoreTot_temp = zeros(numWorkers,1);
randOuterCScoreTot_temp = zeros(numWorkers,1);
skipOuterCScoreTot_temp = zeros(numWorkers,1);
randSkipOuterCScoreTot_temp = zeros(numWorkers,1);
randCnt_temp = zeros(numWorkers,1);
skipRandCnt_temp = zeros(numWorkers,1);
innerPreRandRxnCntTot_temp = zeros(numWorkers,1);
randRxnCntTot_temp = zeros(numWorkers,1);
innerRandRxnCntTot_temp = zeros(numWorkers,1);
innerSkipRxnCntTot_temp = zeros(numWorkers,1);
randSkipRxnCntTot_temp = zeros(numWorkers,1);
innerRandSkipRxnCntTot_temp = zeros(numWorkers,1);
innerForLoopCntTot_temp = zeros(numWorkers,1);
innerLoopCntTot_temp = zeros(numWorkers,1);
preRandCScoreTot_temp = zeros(numWorkers,1);
randCScoreTot_temp = zeros(numWorkers,1);
skipCScoreTot_temp = zeros(numWorkers,1);
randSkipCScoreTot_temp = zeros(numWorkers,1);
exConfFlag_temp = zeros(numWorkers,1);
exTrialFlag_temp = zeros(numWorkers,1);
%possList_temp = zeros(numWorkers,1);

%preRandRxnCnt= zeros(numWorkers,1);
%randRxnCnt = zeros(numWorkers,1);
%skipRxnCnt = zeros(numWorkers,1);
randSkipRxnCnt = zeros(numWorkers,1);
innerPreRandRxnCnt = zeros(numWorkers,1);
innerRandRxnCnt = zeros(numWorkers,1);
innerSkipRxnCnt = zeros(numWorkers,1);
innerRandSkipRxnCnt = zeros(numWorkers,1);
%newMetsCnt= zeros(numWorkers,1); 
randNewMetsCnt = zeros(numWorkers,1); 
skipNewMetsCnt = zeros(numWorkers,1);
randSkipNewMetsCnt = zeros(numWorkers,1);
%newMetsProdCnt = zeros(numWorkers,1);
skipNewMetsProdCnt = zeros(numWorkers,1);
innerNewMetsCnt = zeros(numWorkers,1);
innerRandNewMetsCnt = zeros(numWorkers,1);
innerSkipNewMetsCnt = zeros(numWorkers,1);
innerRandSkipNewMetsCnt = zeros(numWorkers,1);
%outerCScore = zeros(numWorkers,1);
%innerLoopCnt = zeros(numWorkers,1);
%randRxnCnt = zeros(numWorkers,1);
%skipRxnCnt = zeros(numWorkers,1);
randSkipRxnCnt = zeros(numWorkers,1);
innerRandRxnCnt = zeros(numWorkers,1); 
innerSkipRxnCnt = zeros(numWorkers,1);
%newMetsCnt = zeros(numWorkers,1); 
randNewMetsCnt = zeros(numWorkers,1);
skipNewMetsCnt = zeros(numWorkers,1);
randSkipNewMetsCnt = zeros(numWorkers,1);
%newMetsProdCnt = zeros(numWorkers,1);
skipNewMetsProdCnt = zeros(numWorkers,1);
innerNewMetsCnt = zeros(numWorkers,1);
innerRandNewMetsCnt = zeros(numWorkers,1);
innerSkipNewMetsCnt = zeros(numWorkers,1);
innerRandSkipNewMetsCnt = zeros(numWorkers,1);
%outerCScore = zeros(numWorkers,1);
%innerLoopCnt = zeros(numWorkers,1);
%pathLength = zeros(numWorkers,1);



%dataStore_temp = strings{numWorkers,1,50};
% trialInfo_temp = strings(numWorkers,1,25);

%ensure to check met list can be divided into equal groups for parallel processing
%by adding dummy metabolites to ensure a remainder of 0 when seperating by
%number 

rc =rem(length(checkMets),numWorkers);
while rc > 0 
    checkMets = [checkMets;"cpd99999[c]"];
    rc = rc-1;
end

checkDistSize = floor(length(checkMets)/numWorkers);
checkMets_temp = strings(numWorkers,checkDistSize);

for W = 1:numWorkers
    %Divide into seperate pools for parallel processing
    if W == 1
        checkMets_temp(W,:) = checkMets(1:(checkDistSize*W));
    %elseif W == numWorkers
        %checkMets_temp(W,:) = checkMets((checkDistSize*(W-1)):end);
    else
        checkMets_temp(W,:) = checkMets((checkDistSize*(W-1)+1):(checkDistSize*W));
    end
    outerLoopCnt_temp(W,1) = outerLoopCnt;
    exTrialFlag_temp(W,1) = exTrialFlag;
    exConfFlag_temp(W,1) = exConfFlag;
    masterFlag_temp(W,1) = masterFlag;
    largePathFlag_temp(W,1) = largePathFlag;
    viableMets_temp(W,:) = viableMets;
    blockedMets_temp(W,:) = blockedMets;
    
    avgPoss_temp(W,1) = avgPoss;
    trial_temp(W,1) = trial;
    skipCnt_temp(W,1) = skipCnt;
    randSkipCnt_temp(W,1) = randSkipCnt;
    outeridx_temp(W,1) = outeridx;
    maxTrials_temp(W,1) = maxTrials;
    maxLength_temp(W,1) = maxLength;
    innerRandFlag_temp(W,1) = innerRandFlag;
    idx_shift_temp(W,1) = idx_shift;
    preRandRxnCntTot_temp(W,1) = preRandRxnCntTot;
    skipRxnCntTot_temp(W,1) = skipRxnCntTot;
    outerCScoreTot_temp(W,1) = outerCScoreTot;
    randOuterCScoreTot_temp(W,1) = randOuterCScoreTot;
    skipOuterCScoreTot_temp(W,1) = skipOuterCScoreTot;
    randSkipOuterCScoreTot_temp(W,1) = randSkipOuterCScoreTot;
    randCnt_temp(W,1) = randCnt;
    skipRandCnt_temp(W,1) = skipRandCnt;
    innerPreRandRxnCntTot_temp(W,1) = innerPreRandRxnCntTot;
    randRxnCntTot_temp(W,1) = randRxnCntTot;
    innerRandRxnCntTot_temp(W,1) = innerRandRxnCntTot;
    innerSkipRxnCntTot_temp(W,1) = innerSkipRxnCntTot;
    randSkipRxnCntTot_temp(W,1) = randSkipRxnCntTot;
    innerRandSkipRxnCntTot_temp(W,1) = innerRandSkipRxnCntTot;
    innerForLoopCntTot_temp(W,1) = innerForLoopCntTot;
    innerLoopCntTot_temp(W,1) = innerLoopCntTot;
    preRandCScoreTot_temp(W,1) = preRandCScoreTot;
    randCScoreTot_temp(W,1) = randCScoreTot;
    skipCScoreTot_temp(W,1) =  skipCScoreTot;
    randSkipCScoreTot_temp(W,1) = randSkipCScoreTot;
    exConfFlag_temp(W,1) = exConfFlag;
    exTrialFlag_temp(W,1) = exTrialFlag;
    %possList = possList;
%     trialInfo_temp{W,1,25} = trialInfo;
    dataStore_temp = dataStore;
    
    %preRandRxnCnt = preRandRxnCnt;
    %randRxnCnt = randRxnCnt;
    %skipRxnCnt = skipRxnCnt;
    randskipRxnCnt = randSkipRxnCnt;
    innerPreRandRxnCnt = innerPreRandRxnCnt;
    innerrandRxnCnt = innerRandRxnCnt;
    innerskipRxnCnt = innerSkipRxnCnt;
    innerRandSkipRxnCnt = innerRandSkipRxnCnt;
    newMetsCnt = newMetsCnt; 
    randnewMetsCnt = randNewMetsCnt; 
    skipnewMetsCnt = skipNewMetsCnt;
    randSkipnewMetsCnt = randSkipNewMetsCnt;
    %newMetsProdCnt = newMetsProdCnt;
    skipnewMetsProdCnt = skipNewMetsProdCnt;
    innernewMetsCnt = innerNewMetsCnt;
    innerRandnewMetsCnt = innerRandNewMetsCnt;
    innerSkipnewMetsCnt = innerSkipNewMetsCnt;
    innerRandSkipnewMetsCnt = innerRandSkipNewMetsCnt;
    %outerCScore = outerCScore;
    %innerLoopCnt = innerLoopCnt;
    randRxnCnt = randRxnCnt;
    %skipRxnCnt = skipRxnCnt;
    randskipRxnCnt = randSkipRxnCnt;
    innerrandRxnCnt = innerRandRxnCnt; 
    innerskipRxnCnt = innerSkipRxnCnt;
    %newMetsCnt = newMetsCnt; 
    randnewMetsCnt =  randNewMetsCnt;
    skipnewMetsCnt = skipNewMetsCnt;
    randSkipnewMetsCnt = randSkipNewMetsCnt;
    %newMetsProdCnt = newMetsProdCnt;
    skipnewMetsProdCnt = skipNewMetsProdCnt;
    innernewMetsCnt = innerNewMetsCnt;
    innerRandnewMetsCnt = innerRandNewMetsCnt;
    innerSkipnewMetsCnt = innerSkipNewMetsCnt;
    innerRandSkipnewMetsCnt = innerRandSkipNewMetsCnt;
    %outerCScore = outerCScore;
    innerLoopCnt = innerLoopCnt;
    %pathLength = pathLength;
    
    
    
    
    
end
              
parfor W = 1:numWorkers
    possList=[];
    dataStore_temp_edit = dataStore_temp{W,:,:}
    while ~isempty(checkMets_temp(W,:))
        %Assigning value and initializing values
        cscore  = 0;
        %%%%trialInfo{end+1,1} = initMet;
        gpr_match_origins = {};

        gprInfo = {};
        tempGprInfo = {};
        fullRxnIDs = "";
        outerRxnIDs = "";       %Stores rxnIDs for rxns from outer loop
        innerRxnIDs = "";       %Stores rxnIDs for rxns from inner loop
        fullTracers = "";       %Stores all tracer symbols for trial
        outerTracers = "";      %Stores tracer symbols for rxns from outer loop
        innerTracers = "";      %Stores tracer symbols for rxns from inner loop
        fullRxns = "";          %Stores all rxns for trial
        outerRxns = "";         %Stores rxns from outer loop
        innerRxns = "";         %Stores rxns from inner loop
        innerAllRxns = "";      %Stores all rxns that the inner loop looks through, including those that may lead to infeasible loops
        fullRevs = "";          %Stores all reversibilities of trial
        outerRevs = "";         %Stores reversibilities of rxns from outer loop
        innerRevs = "";         %Stores reversibilities of rxns from inner loop
        fullComps = "";         %Stores all compounds of rxns of trial
        outerComps = "";        %Stores compounds of rxns from outer loop
        innerComps = "";        %Stores compounds of rxns from inner loop
        fullTraceNums = "";     %Stores all trace numbers of trial
        outerTraceNums = "";    %Stores trace values of rxns from outer loop
        innerTraceNums = "";    %Stores  trace values of rxns from inner loop
        ruleBreakLvl = 0;       %Stores level of rule breaking for trial
        skipAddActive = 0;      %Stores if skipAdd is active
        randPermActive = 0;     %Stores if randPerm is active
        trialSuccess = 0;       %Stores fail/success/skipped met/failed met of trial
        failConfCheck = 0;      %Stores if trial failed due to confidence score



        %%%Incrementing counter for number of times through outer loop
        outerLoopCnt_temp(W,1) = outerLoopCnt_temp(W,1) + 1;

        %%%Initializing/Resetting tracker variables
        preRandRxnCnt = 0;              %Rxns added w/ out randPerm and skipAdd
        randRxnCnt = 0;                 %Rxns added w/ randPerm
        skipRxnCnt = 0;                 %Rxns added w/ skipAdd
        randSkipRxnCnt = 0;             %Rxns added w/ randPerm and skipAdd
        innerPreRandRxnCnt = 0;         %Rxns added from inner loop w/ out randPerm and skipAdd
        innerRandRxnCnt = 0;            %Rxns added from inner loop w/ randPerm
        innerSkipRxnCnt = 0;            %Rxns added from inner loop w/ skipAdd
        innerRandSkipRxnCnt = 0;        %Rxns added from inner loop w/ randPerm and skipAdd
        randSkipActive = 0;



        %Definition:
        %new  met: a metabolite that was previously blocked that is added to
        %solution

        newMetsCnt = 0;                 %Tracks new mets w/out randPerm or skipAdd active
        randNewMetsCnt = 0;             %Tracks new mets as a result of randPerm
        skipNewMetsCnt = 0;             %Tracks new mets as a result skipAdd
        randSkipNewMetsCnt = 0;         %Tracks new mets as a result of randPerm and skipAdd
        newMetsProdCnt = 0;             %Tracks new mets that are products added w/out skipAdd
        skipNewMetsProdCnt = 0;         %Tracks new mets that are products added as a result of skipAdd
        innerNewMetsCnt = 0;            %Tracks new mets added in inner loop w/out randPerm and skipAdd
        innerRandNewMetsCnt = 0;        %Tracks new mets added in inner loop w/ randPerm
        innerSkipNewMetsCnt = 0;        %Tracks new mets added in inner loop w/ skipAdd
        innerRandSkipNewMetsCnt = 0;    %Tracks new mets added in inner loop w/ randPerm and skipAdd
        outerCScore = 0;                %Tracks conf score contribution of outer loop to conf score of solution
        innerLoopCnt = 0;               %Tracks number of times through inner loop for solution


        pathLength = 0;                 %stores the path length which is equal to the number of downstream reactions so far added to find a solution
        escapeFlag = 0;                 %used to force exit the outer loop
        solvedFlag = 0;                 %tells the loop to exit and also save and display the found solution
        skipFlag = 0;                   %flags when the reason for flux simulation not working is unknown because all required metabolites are presumably viable, exits the loop
        forceStop = 0;                  %this might not do anything anymore 
        traceNum = 1;                   %number used to track where reactions and metabolites lie downstream from the initial metabolite
        checkMets_temp_init = checkMets_temp(W,:)
        initMet = checkMets_temp_init(end);%initial metabolite whose blockage is being solved for
        initMetName = [];
        
        disp(initMet)
        [initMetName] = retrieveMetName(initMet, masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
        if iscell(initMetName)
            initMetName = initMetName{1};
        end
        if initMet == initMetName
            initMetName = [];
        end
%         trialInfo_temp_edit = trialInfo_temp{W,:,:}
%         trialInfo_temp_edit{end+1,1} = initMet;
%         trialInfo_temp{W,:,:} = trialInfo_temp_edit;


        %Storing rulebreaking level applied
        if ~exTrialFlag_temp(W,1) && ~exConfFlag_temp(W,1) && ~masterFlag_temp(W,1) && ~largePathFlag_temp(W,1)
            ruleBreakLvl = 0;
        elseif exTrialFlag_temp(W,1) && ~exConfFlag_temp(W,1) && ~masterFlag_temp(W,1) && ~largePathFlag_temp(W,1)
            ruleBreakLvl = 1;
        elseif exConfFlag_temp(W,1) && ~masterFlag_temp(W,1) && ~largePathFlag_temp(W,1)
            ruleBreakLvl = 2;
        elseif masterFlag_temp(W,1) && ~largePathFlag_temp(W,1)
            ruleBreakLvl = 3;
        elseif largePathFlag_temp(W,1)
            ruleBreakLvl = 4;
        end

        met = initMet; %met to be passed through additional functions, starts as the initial but will change as the solver goes downstream
        tempBlockMets = blockedMets_temp(W,:);
        %these temporary lists provide a way for the true lists to update as
        %new reactions and metabolites are solved for and added. Once a set of
        %reactions is passed as a solution these lists are used to update the
        %master list, otherwise they are discarded
        tempViaMets = viableMets_temp(W,:);
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
        rulelist = ''; %list of characters that tell the function how to format: c = present compound, s = target/blocked compound, n = new reaction, x = terminator

        %!FIX
        if trial_temp(W,1) > stdTrials/2 || trial_temp(W,1) > avgPoss_temp(W,1) %rem(trial,2)==0 %skipAdd tells the tracerLoop to skip checking for existing reactions first, I think I was probably trying to force skip for even values to diversify the potential solutions but this needs be handled better. Regardless, the ability to skip checking existing reactions at the begining of the search is critical for covering all posibilities
            skipAdd=1;
            skipAddActive = 1;
            randSkipActive = 2; %Signifies that just skipAdd is active for trial, may be changed to the value of 3 if randPerm is activated
            %%%Incrementing appropriate values
            if ~innerRandFlag_temp(W,1)
                skipCnt_temp(W,1) = skipCnt_temp(W,1) + 1;
            else
                randSkipCnt_temp(W,1) = randSkipCnt_temp(W,1) + 1;
            end
        else
            skipAdd=0;
        end
        skipAdd = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract initial reactions, iterate through until end until suitible 
        %solution is found if exhausted all iteretions before then, randomize
        %index for inner indicies
        [addRxns] = tracerLoop2(met,traceNum,tempViaMets,tempExcludeRxns,presentDB,superDB,orthoDB,masterDB,largePathFlag_temp(W,1),masterFlag_temp(W,1),netRxnsMax,skipAdd,presRxns);
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
            elseif ~largePathFlag_temp(W,1) && largePathTemp
                warning("FATALERROR");
                largePathFlag_temp(W,1) = 1;
                outeridx_temp(W,1) = 0;
                maxTrials_temp(W,1)=largePathTrials;
            elseif largePathFlag_temp(W,1) 
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

            if~largePathFlag_temp(W,1) && largePathTemp
                fprintf("No downstream reactions found under standard requirements, reattempting with less stringent parameters\n");
                largePathFlag_temp(W,1) = 1;
                maxLength_temp(W,1)=largeLenMax;
                outeridx_temp(W,1) = 0;
            elseif largePathFlag_temp(W,1)
                fprintf("No downstream reactions found under easement, temporarily disabling metabolite from the blocked list\nAdding metabolite to seperate output list, consider manual validation of associated reaction(s)\n");
                escapeFlag = 1;
            else
                fprintf("No downstream reactions found under standard requirements,passing\n");
                escapeFlag=1;
            end
        else
            nextLen=length(addRxns(:,1)); %number of qualifying reactions found
            possib=nextLen; %number of possibilities to test
            if nextLen>outeridx_temp(W,1) %iterates through all possible solutions 
                outeridx_temp(W,1)=outeridx_temp(W,1)+1;
            else
                %!FIX
                outeridx_temp(W,1) = 1; %after trying all potential solutions for the first reaction, the selected first solution is randomized -> I do not like this approach an alternative should be found here
                innerRandFlag_temp(W,1) = 1;
                idx_shift_temp(W,1) = idx_shift_temp(W,1)-1;
                randPermActive = 1; %stores if randPerm is active for the trialInfo output file
                %Stores if randPerm and skipAdd are active (assigns value of 3)
                %or just randPerm is active (assigns value of 1)
                if ~skipAdd
                    randSkipActive = 3;
                else
                    randSkipActive = 1;
                end
            end

            %after a reaction is selected, its information is extracted into
            %the followingl lists which are used in the inner loop where the
            %next reaction will be selected using a similar algorithm
            rxnID = addRxns{outeridx_temp(W,1),1}; 
            react = transpose(addRxns{outeridx_temp(W,1),2});
            prodt = transpose(addRxns{outeridx_temp(W,1),3});
            rxn = addRxns{outeridx_temp(W,1),4};
            rev = addRxns{outeridx_temp(W,1),5};
            trace = addRxns{outeridx_temp(W,1),6};
            %Saving gene information, so GPR match can be dtermined if a
            %suitable solution is found
            trNumber = addRxns{outeridx_temp(W,1),7}; 
    %         if length(addRxns(1,:)) == 8
    %             gpr = addRxns{outeridx,8};
    %         else
    %             gpr = [];
    %         end
            %Saving gene origins if reaction has a GPR match
            if trace == "(+)"
                gpr_match_origins{end+1, 1} = addRxns{outeridx_temp(W,1), 8};
                gpr_match_origins{end, 2} = addRxns{outeridx_temp(W,1), 9};
                gpr_match_origins{end, 3} = addRxns{outeridx_temp(W,1), 10};
            end
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
                    outerComps = outerComps + "(*" + prodt(i) + "*)";
                elseif ismember(prodt(i),tempViaMets)
                    [metName] = retrieveMetName(prodt(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                    clist = [clist;metName];
                    outerComps = outerComps + "(" + prodt(i) + ")";
                    rulelist(end+1) = 's';
                    ccnt = ccnt+1;
                else
                    [metName] = retrieveMetName(prodt(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                    clist = [clist;metName+"*"];
                    outerComps = outerComps + "(" + prodt(i) + "*)";
                    rulelist(end+1) = 'c';
                    scnt=scnt+1;
                    tempViaMets=[tempViaMets;prodt(i)];
                    %%%Incrementing appropriate values
                    if ~skipAdd
                        newMetsCnt = newMetsCnt + 1;
                        newMetsProdCnt = newMetsProdCnt + 1;
                    else
                        skipNewMetsCnt = skipNewMetsCnt + 1;
                        skipNewMetsProdCnt = skipNewMetsProdCnt + 1;
                    end
                end
            end
            rulelist(end+1) = 'n';
            for i=1:length(react)
                currentMets=[currentMets;react(i)];
                if ismember(react(i),viableMets_temp(W,:))
                    [metName] = retrieveMetName(react(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                    clist = [clist;metName];
                    metIdx = find(contains(masterMetIDs, react(i)));
                    outerComps = outerComps + "(" + react(i) + ")";
                    rulelist(end+1) = 's';
                    ccnt = ccnt+1;
                else
                    [metName] = retrieveMetName(react(i), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                    clist = [clist;metName+"*"];
                    outerComps = outerComps + "(" + react(i) + "*)";
                    rulelist(end+1) = 'c';
                    loopMets = [loopMets;react(i)];
                    scnt=scnt+1;
                    tempViaMets=[tempViaMets;react(i)];
                    %%%Incrementing appropriate values
                    if ~skipAdd
                        newMetsCnt = newMetsCnt + 1;
                    else
                        skipNewMetsCnt = skipNewMetsCnt + 1;
                    end
                end
            end
            rulelist(end+1) = 'n';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update lists of reaction information for visual translation,
            %store information not already contained in the model in temporary
            %lists to be used if a solution is found
            tlist = [tlist;trace];
    %         gprInfo{end+1} = gpr;
            %outerTracers = outerTracers + "[" + trace + "]";
            outerTracers = outerTracers + "[" + trace + "]";
            rklist = [rklist;trNumber];
            outerTraceNums = outerTraceNums + "(" + num2str(trNumber) + ")";
            revlist = [revlist;rev];  
            outerRevs = outerRevs + "(" + num2str(rev) + ")";
            rxlist = [rxlist;rxn];
            outerRxns = outerRxns + "(" + rxn + ")";

            tempExcludeRxns = [tempExcludeRxns;rxnID];   


            if ~any(strcmp(presID,rxnID))
                tempAddID = [tempAddID;rxnID];
                tempAddRxn = [tempAddRxn;rxn];
                tempAddRev = [tempAddRev;rev];
                tempTrace = [tempTrace;trace];
    %             tempGprInfo{end+1} = gpr;
                rlist = [rlist;rxnID+"*"];
                outerRxnIDs = outerRxnIDs + "(" + rxnID +  "*)";
                rcnt=rcnt+1;
            else
                rlist = [rlist;rxnID];

                outerRxnIDs = outerRxnIDs + "(" + rxnID + ")";
            end
            %%%Incrementing appropriate counters
            if ~skipAdd
                preRandRxnCnt = preRandRxnCnt + 1;
                preRandRxnCntTot_temp(W,1) = preRandRxnCntTot_temp(W,1) + 1;
            elseif skipAdd
                skipRxnCnt = skipRxnCnt + 1;
                skipRxnCntTot_temp(W,1) = skipRxnCntTot_temp(W,1) + 1;
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

            %%%Incrementing appropriate values -- confidence score contribution
            %%%of outer loop
    %         totgene = 0;
    %         totX = 0;
    %         totT = 0;
    %         tot0 = 0;
    %         tot1 = 0;
    %         tot2 = 0;
    %         tot3 = 0;
    %         tot4 = 0;
    %         tot5 = 0;
    %         tot6 = 0;
            for i = 1:length(tlist)
                sym = string(extractBetween(tlist(i),"(",")"));
                if sym == "+"
                    outerCScore = outerCScore + originData(1).confScore; 
                elseif sym == "X"
                    outerCScore = outerCScore + originData(2).confScore; 
                elseif sym == "T"
                    outerCScore = outerCScore + originData(3).confScore; 
                elseif sym == "0"
                    outerCScore = outerCScore + originData(4).confScore; 
    %             elseif sym == "1"
    %                 tot1 = tot1+1;
    %             elseif sym == "2"
    %                 tot2 = tot2+1;
    %             elseif sym == "3"
    %                 tot3 = tot3+1;
    %             elseif sym == "4"
    %                 tot4 = tot4+1;
    %             elseif sym == "5"
    %                 tot5 = tot5+1;
    %             elseif sym == "6"
    %                 tot6 = tot6+1;
                elseif str2double(sym)+4 <= length(originData)
                    try
                        outerCScore = outerCScore + originData(str2double(sym)+4).confScore;
                    catch
                        error("INVALID TRACER SYMBOL\nRank is not a proper index");
                    end
                elseif str2double(sym)+4 > length(originData)
                    error("INVALID TRACER SYMBOL\nRank exceeds those given.");
                else
                    error("INVALID TRACER SYMBOL\nSymbols must be formated such that the traceNum is directly between parenthesis (ex. '<(#)>' ) \n\n");
                end
            end
            %Calculation of Confidence Score
    %         cscore = scoreCalc(totgene,totX,totT,tot0,tot1,tot2,tot3,tot4,tot5,tot6,rcnt,ccnt);
    %         outerCScore = cscore;
            if ~innerRandFlag_temp(W,1) && ~skipAdd
                outerCScoreTot_temp(W,1) = outerCScoreTot_temp(W,1) + outerCscore;
            elseif ~skipAdd
                randOuterCScoreTot_temp(W,1) = randOuterCScoreTot_temp(W,1) + outerCScore;
            elseif ~innerRandFlag_temp(W,1)
                skipOuterCScoreTot_temp(W,1) = skipOuterCScoreTot_temp(W,1) + outerCScore;
            else
                randSkipOuterCScoreTot_temp(W,1) = randSkipOuterCScoreTot_temp(W,1) + outerCScore;
            end

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
            while ~isempty(loopMets) && ~forceStop && pathLength<=maxLength_temp(W,1)
                lastLen=nextLen; 
                traceNum = traceNum+1;
                newLoopMets = string([]);
                for i = 1:length(loopMets)
                    met = loopMets(i);
                    if met ~= "X"
                        rulelist(end) = 'n';
                    end
                    [addRxns] = tracerLoop2(met,traceNum,tempViaMets,tempExcludeRxns,presentDB,superDB,orthoDB,masterDB,largePathFlag,masterFlag,netRxnsMax,skipAdd,presRxns);
                    %Randomize Index after nonrandomized attempts fail            
                    if iscell(addRxns) && ~isempty(addRxns) && pathLength<maxLength_temp(W,1)
                        pcheck = char(addRxns{1,6}); %checks the origin of the first/best qualifying reaction from tracerLoop
                        inTrialList=[];
                        for j=1:length(addRxns(:,1))
                            inTrialList=[inTrialList;j];
                        end
                        %!FIX - do not want function to be dependent on randomization
                        %??? Should be and? 
                        if innerRandFlag_temp(W,1) %%%&& contains(pcheck,'0')%%% %if the randomizer flag from before is active (meaning no solution was found for all of the outerloop designated starting point reaction) AND the first/best qualifying reations for the next reaction is a reaction already in the model,
                            %idx_shift = idx_shift - 1;
                            inTrialList = circshift(inTrialList, idx_shift_temp(W,1));
                            %%%inTrialList=inTrialList(randperm(length(inTrialList))); %then the list of inner loop reactions to test will be randomized - this attempts to allow for more comb inations but also may be a massive time waster so optimization here should be looked into
                            %%%Incrementing appropriate counter
                            if ~skipAdd
                                randCnt_temp(W,1) = randCnt_temp(W,1) + 1;
                            else
                                skipRandCnt_temp(W,1) = skipRandCnt_temp(W,1) + 1; 
                            end

                        end

                        chidx=1;
                        passFlag=0;
                        nextLen=length(inTrialList);
                        possib= nextLen*lastLen; %total number of possibilities
    %                     warning('stop here')
                        while chidx<=nextLen && ~passFlag
                                innertrial=inTrialList(chidx);
                                chidx=chidx+1;
                                rxn = addRxns{innertrial,4};
                                innerAllRxns = [innerAllRxns+"(" + rxn + ")"];
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
                                %If reaction has a GPR match, saving origins
                                %and corresponding genes
                                if trace == "(+)"
                                    gpr_match_origins{end+1, 1} = addRxns{innertrial, 8};
                                    gpr_match_origins{end, 2} = addRxns{innertrial, 9};
                                    gpr_match_origins{end, 3} = addRxns{innertrial, 10};
                                end
                                trNumber = addRxns{innertrial,7};
    %                             if length(addRxns(1,:)) == 8
    %                                 gpr = addRxns{innertrial,8};
    %                             else
    %                                 gpr = [];
    %                             end
                                passFlag=1;
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Update lists for visual translation
                        if passFlag
                            for j=1:length(react)
                                currentMets=[currentMets;react(j)];                      
                                if  ismember(react(j),viableMets_temp(W,:))
                                    [metName] = retrieveMetName(react(j), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                                    clist = [clist;metName];
                                    innerComps = innerComps + "(" + react(j) + ")";
                                    rulelist(end+1) = 's';
                                    scnt = scnt+1;
                                else
                                    [metName] = retrieveMetName(react(j), masterMets, masterMetsOld, masterMetIDs, masterMetIDsOld);
                                    clist = [clist;metName+"*"];
                                    metIdx = find(contains(masterMetIDs, react(j)));
                                    innerComps = innerComps + "(" + react(j) + ")";
                                    newLoopMets = [newLoopMets;react(j)];
                                    rulelist(end+1) = 'c';
                                    ccnt = ccnt+1;
                                    tempViaMets=[tempViaMets;react(j)];
                                    %%%Incrementing appropriate values
                                    if ~innerRandFlag_temp(W,1) && ~skipAdd
                                        newMetsCnt = newMetsCnt + 1; 
                                        innerNewMetsCnt = innerNewMetsCnt + 1;
                                    elseif ~skipAdd
                                        randNewMetsCnt = randNewMetsCnt + 1;
                                        innerRandNewMetsCnt = innerRandNewMetsCnt + 1;
                                    elseif ~innerRandFlag_temp(W,1)
                                        skipNewMetsCnt = skipNewMetsCnt + 1; 
                                        innerSkipNewMetsCnt = innerSkipNewMetsCnt + 1; 
                                    else
                                        randSkipNewMetsCnt = randSkipNewMetsCnt + 1;
                                        innerRandSkipNewMetsCnt = innerRandSkipNewMetsCnt + 1;
                                    end
                                end
                            end
                            clist = [clist;"X"];
                            tlist = [tlist;trace];
        %                     gprInfo{end+1} = gpr;
                            innerTracers = innerTracers + "[" + trace + "]";
                            revlist = [revlist;rev];
                            innerRevs = innerRevs + "(" + num2str(rev) + ")";
                            rxlist = [rxlist;rxn];
                            innerRxns = innerRxns + "(" + rxn + ")";
                            innerTraceNums = innerTraceNums + "(" + num2str(trNumber) + ")";
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            tempExcludeRxns = [tempExcludeRxns;rxnID];
                            if ~any(strcmp(presID,rxnID))
                                tempAddID = [tempAddID;rxnID];
                                tempAddRxn = [tempAddRxn;rxn];
                                tempAddRev = [tempAddRev;rev];
                                tempTrace = [tempTrace;trace];
        %                         tempGprInfo{end+1} = gpr;

                                rlist = [rlist;rxnID+"*"];
                                innerRxnIDs = innerRxnIDs + "(" + rxnID + "*)";
                                rcnt = rcnt+1;
                            else
                                rlist = [rlist;rxnID];

                                innerRxnIDs = innerRxnIDs + "(" + rxnID + ")";
                            end
                            %%%Incrementing appropriate counters
                            if ~innerRandFlag_temp(W,1) && ~skipAdd
                                preRandRxnCnt = preRandRxnCnt + 1;
                                preRandRxnCntTot_temp(W,1) = preRandRxnCntTot_temp(W,1) + 1;
                                innerPreRandRxnCnt = innerPreRandRxnCnt + 1;
                                innerPreRandRxnCntTot_temp(W,1) = innerPreRandRxnCntTot_temp(W,1) + 1;
                            elseif ~skipAdd
                                randRxnCnt = randRxnCnt + 1;
                                randRxnCntTot_temp(W,1) = randRxnCntTot_temp(W,1) + 1;
                                innerRandRxnCnt = innerRandRxnCnt + 1;
                                innerRandRxnCntTot_temp(W,1) = innerRandRxnCntTot_temp(W,1) + 1;
                            elseif ~innerRandFlag_temp(W,1)
                                skipRxnCnt = skipRxnCnt + 1;
                                skipRxnCntTot_temp(W,1) = skipRxnCntTot_temp(W,1) + 1;
                                innerSkipRxnCnt = innerSkipRxnCnt + 1;
                                innerSkipRxnCntTot_temp(W,1) = innerSkipRxnCntTot_temp(W,1) + 1;
                            else
                                randSkipRxnCnt = randSkipRxnCnt + 1;
                                randSkipRxnCntTot_temp(W,1) = randSkipRxnCntTot_temp(W,1) + 1;
                                innerRandSkipRxnCnt = innerRandSkipRxnCnt + 1;
                                innerRandSkipRxnCntTot_temp(W,1) = innerRandSkipRxnCntTot_temp(W,1) + 1;
                            end
                        end
                    else
                        forceStop = 1;
                        rulelist(end+1) = 'f';
                    end
                innerForLoopCntTot_temp(W,1) = innerForLoopCntTot_temp(W,1) + 1;    
                end
                checkMet_temp_flip = checkMets_temp(W,:)
                checkMets_temp_flip = flipud(checkMets_temp_flip);
                checkMets_temp_flip = flipud([checkMets_temp_flip;"X"]);
                
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
                innerLoopCnt = innerLoopCnt + 1;
            end
            innerLoopCntTot_temp(W,1) = innerLoopCntTot_temp(W,1) + innerLoopCnt;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %If a pathway is sucessfully solved to completion, calculate its 
            %confidence score. If acceptable solution is found, escape the loop,
            %reset starting parameters,commit changes to blocked/viable lists 
            %and add rxns + corresponding IDs to output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Summation of tracer values
            if rulelist(end) == 'x'%         totgene = 0;
        %         totX = 0;
        %         totT = 0;
        %         tot0 = 0;
        %         tot1 = 0;
        %         tot2 = 0;
        %         tot3 = 0;
        %         tot4 = 0;
        %         tot5 = 0;
        %         tot6 = 0;
                for i = 1:length(tlist)
                    sym = string(extractBetween(tlist(i),"(",")"));
                    if sym == "+"
                        cscore = cscore + originData(1).confScore; 
                    elseif sym == "X"
                        cscore = cscore + originData(2).confScore; 
                    elseif sym == "T"
                        cscore = cscore + originData(3).confScore; 
                    elseif sym == "0"
                        cscore = cscore + originData(4).confScore; 
        %             elseif sym == "1"
        %                 tot1 = tot1+1;
        %             elseif sym == "2"
        %                 tot2 = tot2+1;
        %             elseif sym == "3"
        %                 tot3 = tot3+1;
        %             elseif sym == "4"
        %                 tot4 = tot4+1;
        %             elseif sym == "5"
        %                 tot5 = tot5+1;
        %             elseif sym == "6"
        %                 tot6 = tot6+1;
                    elseif str2double(sym)+4 <= length(originData)
                        try
                            outerCScore = outerCScore + originData(str2double(sym)+4).confScore;
                        catch
                            error("INVALID TRACER SYMBOL\nRank is not a proper index");
                        end
                    elseif str2double(sym)+4 > length(originData)
                        error("INVALID TRACER SYMBOL\nRank exceeds those given.");
                    else
                        error("INVALID TRACER SYMBOL\nSymbols must be formated such that the traceNum is directly between parenthesis (ex. '<(#)>' ) \n\n");
                    end
                end
                %Calculation of Confidence Score
    %             warning('Look here')
    %             cscore = scoreCalc(totgene,totX,totT,tot0,tot1,tot2,tot3,tot4,tot5,tot6,rcnt,ccnt);

                %%%Incrementing appropriate variables
                if ~innerRandFlag_temp(W,1) && ~skipAdd
                    preRandCScoreTot_temp(W,1) = preRandCScoreTot_temp(W,1) + cscore;
                elseif ~skipAdd
                    randCScoreTot_temp(W,1) = randCScoreTot_temp(W,1) + cscore;
                elseif ~innerRandFlag_temp(W,1)
                    skipCScoreTot_temp(W,1) = skipCScoreTot_temp(W,1) + cscore;
                else
                    randSkipCScoreTot_temp(W,1) = randSkipCScoreTot_temp(W,1) + cscore;
                end

    %             if exTrialFlag && ~easementFlag
    %                 cscore = cscore + extraPenalty(cscore);
    %             elseif easementFlag
    %                 cscore = cscore + easementPenalty(cscore);
    %             end

                %confidence scores that do not pass have the option of being
                %resolved under different parameters specified in the user
                %config at the begining of the script
                if (cscore <= reqConf && ~exConfFlag_temp(W,1) || cscore <=boostConf && exConfFlag_temp(W,1) && ~masterFlag || cscore<=mastConf && masterFlag || cscore<=largePathConf && largePathFlag) && (~all(contains(rlist,presRxns)) || skipAdd)
                    if ~exTrialFlag_temp(W,1) && ~exConfFlag_temp(W,1) && ~masterFlag && ~largePathFlag 
                        fprintf(outFile,"\rRESULTS FOR METABOLITE %s\t%s\n\nSolved under standard parameters\n\n", initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved under standard parameters\n\n",initMet,initMetName);
                    elseif exTrialFlag_temp(W,1) && ~exConfFlag_temp(W,1)
                        fprintf(outFile, "\rRESULTS FOR METABOLITE %s\t%s\n\nSolved given additional permutation\n\n", initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved given additional permutation\n\n", initMet,initMetName);
                    elseif exConfFlag_temp(W,1) && ~masterFlag
                        fprintf(outFile, "\rRESULTS FOR METABOLITE %s\t%s\n\nSolved for increased confidence threshold\n\n",initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved for increased confidence threshold\n\n",initMet,initMetName);
                    elseif masterFlag && ~largePathFlag
                        fprintf(outFile,"\rRESULTS FOR METABOLITE %s\t%s\n\nSolved using Master model\n\n",initMet,initMetName);
                        fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nSolved using Master model\n\n",initMet,initMetName);
                    elseif largePathFlag

                    else
                        %error("\rRESULTS FOR METABOLITE "+initMet+"\n\rSolved under unreasonable parameters\n\r");
                    end
                    if ~all(contains(rlist,presRxns)) || skipAdd
                        escapeFlag = 1;
                        solvedFlag = 1;
                    end
                elseif ~all(contains(rlist,presRxns))
                    failConfCheck = 1;
                end
            end
            %Advance to next trial if failed to produce an acceptable
            %pathway. If max trial is hit,either change parameters or pass
            %metabolite
            if ~escapeFlag
                possList=[possList,possib];
                trial_temp(W,1) = trial_temp(W,1)+1;
                avgPoss_temp(W,1)=mean(possList);
                if mod(trial_temp(W,1),10) == 0
                    fprintf(trial_temp(W,1) +"\n");
    %                  avgPoss=mean(possList);
                    fprintf("\naverage # of Possibilities: %s\n",string(avgPoss_temp(W,1)))
                    %this uses the total number of calculated possibilities
                    %(averaged across all trials) to determine when to force
                    %quit. 2 times the calculated possibilities are used since
                    %randomization is used so this should be updated when that
                    %is fixed. 
                    if trial_temp(W,1)>2*avgPoss_temp(W,1)
                        trial_temp(W,1)=999999;
                        exTrialFlag_temp(W,1) = 1;
                    end             
                elseif skipAdd && trial_temp(W,1)>2*avgPoss_temp(W,1)
                    fprintf(trial_temp(W,1) +"\n");
    %                  avgPoss=mean(possList);
                    fprintf("\naverage # of Possibilities: %s\n",string(avgPoss_temp(W,1)))
                    trial_temp(W,1)=999999;
                    exTrialFlag_temp(W,1) = 1;
                end
                if trial_temp(W,1)>=maxTrials_temp(W,1) 
                    fprintf("\rRESULTS FOR METABOLITE %s\t%s\n\nNo Solutions could be determined in the alotted number of trials,\n",initMet,initMetName);
                    if ~allowExTrials && ~allowExConf && ~allowMaster && ~largePathTemp
                        fprintf("Standard (%d) was insufficient for a complete solution. Consider enabling rulebreaking and/or manual validation of associated reaction(s)\n\n",stdTrials);
                        escapeFlag=1;
                        idx_shift_temp(W,1) = 0;
                    elseif allowExTrials && ~exTrialFlag_temp(W,1)
                        fprintf("Standard (%d) was insufficient for a complete solution.\nReattempting with extra trials\n\rTrial:\n",stdTrials);
                        exTrialFlag_temp(W,1) = 1;
                        maxTrials_temp(W,1) = extraTrials;
                        trial_temp(W,1) = 1;
                        outeridx_temp(W,1)=0;
                        idx_shift_temp(W,1) = 0;
                    elseif allowExConf && ~exConfFlag_temp(W,1)
                        if trial_temp(W,1) ~= 999999
                            fprintf("Increasing max number of trials was still insufficient for a complete solution.\nReattempting with less permissible confidence\n\rTrial:\n");
                        else
                            fprintf("Standard (%d) was insufficient for a complete solution.\nExtra trials will also fail.\nReattempting with less permissible confidence\n\rTrial:\n",stdTrials);
                        end
                        exConfFlag_temp(W,1) = 1;
                        trial_temp(W,1) = 1;
                        outeridx_temp(W,1)=0;    
                        idx_shift_temp(W,1) = 0;
                    elseif allowMaster && ~masterFlag
                        fprintf("Increasing permissible confidence was still insufficient for a complete solution.\nReattempting with use of master database\n\rTrial:\n");
                        masterFlag = 1;
                        trial_temp(W,1) = 1;
                        outeridx_temp(W,1)=0;
                        idx_shift_temp(W,1) = 0;
                    elseif largePathTemp && ~largePathFlag
                        fprintf("Utilization of master list was still insufficient for a complete solution.\nReattempting for larger and/or more open pathways\n\rTrial:\n");
                        outeridx_temp(W,1)=0;
                        trial_temp(W,1)=1;
                        largePathFlag=1;
                        maxTrials_temp(W,1) = largePathTrials;
                        idx_shift_temp(W,1) = 0;
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
           outeridx_temp(W,1) = 0;
           innerRandFlag_temp(W,1) = 0;
           trialStore = trial_temp(W,1);
           trial_temp(W,1) = 1;
           idx_shift_temp(W,1) = 0; 
           if solvedFlag && ~all(contains(rlist,presRxns))
               viableMets_temp(W,:) = tempViaMets;   
               trialSuccess = 1; 

               %%%Adding results as init met is solved
               %Assigns value based on level of rulebreaking allowed

               if ~exTrialFlag_temp(W,1) && ~exConfFlag_temp(W,1) && ~masterFlag && ~largePathFlag %Standard parameters
                    dataStore_temp_edit{end+1, 1} = initMet;
                    dataStore_temp_edit{end, 2} = 0;
                elseif exTrialFlag_temp(W,1) && ~exConfFlag_temp(W,1) %with extra trials
                    dataStore_temp_edit{end+1, 1} = initMet;
                    dataStore_temp_edit{end, 2} = 1;
                elseif exConfFlag_temp(W,1) && ~masterFlag  %with extra conf
                    dataStore_temp_edit{end+1, 1} = initMet;
                    dataStore_temp_edit{end, 2} = 2;
                elseif masterFlag && ~largePathFlag %using masterDB
                    dataStore_temp_edit{end+1, 1} = initMet;
                    dataStore_temp_edit{end, 2} = 3;
                elseif largePathFlag  %with largePath
                    dataStore_temp_edit{end+1, 1} = initMet;
                    dataStore_temp_edit{end, 2} = 4;
                else
                    %error("\rRESULTS FOR METABOLITE "+initMet+"\n\rSolved under unreasonable parameters\n\r");
                    dataStore_temp_edit{end+1, 1} = initMet;
                    dataStore_temp_edit{end, 2} = 5;
                end

                dataStore_temp_edit{end, 3} = trialStore;
                dataStore_temp_edit{end, 4} = pathLength;
                dataStore_temp_edit{end, 5} = randCnt_temp(W,1);
                dataStore_temp_edit{end, 6} = skipRandCnt_temp(W,1);
                dataStore_temp_edit{end, 7} = skipCnt_temp(W,1);
                dataStore_temp_edit{end, 8} = randSkipCnt_temp(W,1);
                dataStore_temp_edit{end, 9} = preRandRxnCnt;
                dataStore_temp_edit{end, 10} = randRxnCnt;
                dataStore_temp_edit{end, 11} = skipRxnCnt;
                dataStore_temp_edit{end, 12} = randSkipRxnCnt;
                dataStore_temp_edit{end, 13} = preRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 14} = randRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 15} = skipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 16} = randSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 17} = innerPreRandRxnCnt;
                dataStore_temp_edit{end, 18} = innerRandRxnCnt;
                dataStore_temp_edit{end, 19} = innerSkipRxnCnt;
                dataStore_temp_edit{end, 20} = innerRandSkipRxnCnt;
                dataStore_temp_edit{end, 21} = innerPreRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 22} = innerRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 23} = innerSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 24} = innerRandSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 25} = newMetsCnt;
                dataStore_temp_edit{end, 26} = randNewMetsCnt;
                dataStore_temp_edit{end, 27} = skipNewMetsCnt;
                dataStore_temp_edit{end, 28} = randSkipNewMetsCnt;
                dataStore_temp_edit{end, 29} = newMetsProdCnt;
                dataStore_temp_edit{end, 30} = skipNewMetsProdCnt;
                dataStore_temp_edit{end, 31} = innerNewMetsCnt;
                dataStore_temp_edit{end, 32} = innerRandNewMetsCnt;
                dataStore_temp_edit{end, 33} = innerSkipNewMetsCnt;
                dataStore_temp_edit{end, 34} = innerRandSkipNewMetsCnt;
                dataStore_temp_edit{end, 35} = cscore;
                dataStore_temp_edit{end, 36} = outerCScore;
                dataStore_temp_edit{end, 37} = outerCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 38} = randOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 39} = skipOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 40} = randSkipOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 41} = preRandCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 42} = randCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 43} = skipCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 44} = randSkipCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 45} = innerLoopCnt;
                dataStore_temp_edit{end, 46} = innerLoopCntTot_temp(W,1);
                dataStore_temp_edit{end, 47} = innerForLoopCntTot_temp(W,1);
                dataStore_temp_edit{end, 48} = outerLoopCnt_temp(W,1);
                dataStore_temp_edit{end, 49} = randSkipActive;
                dataStore_temp_edit{end, 50} = trialSuccess;

                fprintf(dataFile, '%s', initMet);
                fprintf(dataFile, '\t%d', dataStore_temp_edit{end, 2});
                fprintf(dataFile, '\t%d', trialStore);
                fprintf(dataFile, '\t%d', pathLength);
                fprintf(dataFile, '\t%d', randCnt_temp(W,1));
                fprintf(dataFile, '\t%d',skipRandCnt_temp(W,1));
                fprintf(dataFile, '\t%d', skipCnt_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipCnt_temp(W,1));
                fprintf(dataFile, '\t%d', preRandRxnCnt);
                fprintf(dataFile, '\t%d', randRxnCnt);
                fprintf(dataFile, '\t%d', skipRxnCnt);
                fprintf(dataFile, '\t%d', randSkipRxnCnt);
                fprintf(dataFile, '\t%d', preRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', randRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerPreRandRxnCnt);
                fprintf(dataFile, '\t%d', innerRandRxnCnt);
                fprintf(dataFile, '\t%d', innerSkipRxnCnt);
                fprintf(dataFile, '\t%d', innerRandSkipRxnCnt);
                fprintf(dataFile, '\t%d', innerPreRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerRandSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', newMetsCnt);
                fprintf(dataFile, '\t%d', randNewMetsCnt);
                fprintf(dataFile, '\t%d', skipNewMetsCnt);
                fprintf(dataFile, '\t%d', randSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', newMetsProdCnt);
                fprintf(dataFile, '\t%d', skipNewMetsProdCnt);
                fprintf(dataFile, '\t%d', innerNewMetsCnt);
                fprintf(dataFile, '\t%d', innerRandNewMetsCnt);
                fprintf(dataFile, '\t%d', innerSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', innerRandSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', cscore);
                fprintf(dataFile, '\t%d', outerCScore);
                fprintf(dataFile, '\t%d', outerCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', preRandCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerLoopCnt);
                fprintf(dataFile, '\t%d', innerLoopCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerForLoopCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', outerLoopCnt_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipActive);
                fprintf(dataFile, '\t%d\n', trialSuccess);



                randCnt_temp(W,1) = 0;
                skipRandCnt_temp(W,1) = 0;
                skipCnt_temp(W,1) = 0;
                randSkipCnt_temp(W,1) = 0;
                preRandRxnCnt= 0;
                randRxnCnt = 0;
                skipRxnCnt = 0;
                randSkipRxnCnt = 0;
                preRandRxnCntTot_temp(W,1) = 0;
                randRxnCntTot_temp(W,1)= 0;
                skipRxnCntTot_temp(W,1) = 0; 
                randSkipRxnCntTot_temp(W,1) = 0;
                innerPreRandRxnCnt = 0;
                innerRandRxnCnt = 0; 
                innerSkipRxnCnt = 0;
                innerRandSkipRxnCnt = 0;
                innerPreRandRxnCntTot_temp(W,1) = 0;
                innerRandRxnCntTot_temp(W,1) = 0;
                innerSkipRxnCntTot_temp(W,1) = 0;
                innerRandSkipRxnCntTot_temp(W,1) = 0;
                newMetsCnt= 0; 
                randNewMetsCnt = 0; 
                skipNewMetsCnt = 0;
                randSkipNewMetsCnt = 0;
                newMetsProdCnt = 0;
                skipNewMetsProdCnt = 0;
                innerNewMetsCnt = 0;
                innerRandNewMetsCnt = 0;
                innerSkipNewMetsCnt = 0;
                innerRandSkipNewMetsCnt = 0;
    %             cscore = 0;
                outerCScore = 0;
                outerCScoreTot_temp(W,1) = 0;
                randOuterCScoreTot_temp(W,1) = 0;
                skipOuterCScoreTot_temp(W,1) = 0;
                randSkipOuterCScoreTot_temp(W,1) = 0;
                preRandCScoreTot_temp(W,1) = 0;
                randCScoreTot_temp(W,1) = 0;
                skipCScoreTot_temp(W,1) = 0;
                randSkipCScoreTot_temp(W,1) = 0;
                innerLoopCnt = 0;
                innerLoopCntTot_temp(W,1) = 0;
                innerForLoopCntTot_temp(W,1) = 0;
                outerLoopCnt_temp(W,1) = 0;

                %Determining if any reactions have GPR matches
    %             tempGprIdx = 1; %Tracks placement within tempGprInfo list
    %             for m = 1:length(gprInfo)
    %                 if tlist(m) ~= "(+)"
    %                     if ~isempty(gprInfo{m})
    %                         gpr = gprInfo{m};
    %                         [geneInfo] = checkGPRmatch(gpr, orthoDB);
    %                         replaceIdxList = find(strcmp(gpr, gprInfo));
    %                         replaceIdxTemp = find(strcmp(gpr, tempGprInfo));
    %                         for ii = 1:length(replaceIdxList)
    %                             replaceIdx = replaceIdxList(ii);
    %                             gpr_match_origins(replaceIdx,1:3) = geneInfo;
    %                             if ~isempty(geneInfo(1))
    %                                 tlist(replaceIdx) = "(+)";
    %                             end
    %                         end
    %                         for ii = 1:length(replaceIdxTemp)
    %                             replaceIdx = replaceIdxTemp(ii);
    %                             if ~isempty(geneInfo(1))
    %                                 tempTrace(replaceIdx) = "(+)";
    %                             end
    %                         end
    %                     end
    %                 end
    %             end

               if ~masterFlag && ~largePathFlag
                   blockedMets_temp(W,:) = tempBlockMets;
                   checkMets_temp_removeTempBlock = checkMets_temp(W,:)
                   checkMets_temp_removeTempBlock(~contains(checkMets_temp_removeTempBlock,tempBlockMets)) = [];
                   checkMets_temp(W,:) = checkMets_temp_removeTempBlock
                   if exConfFlag_temp(W,1)
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
                   checkMets_temp_removeTempBlock = checkMets_temp(W,:)
                   checkMets_temp(W,:) = checkMets_temp_removeTempBlock(1:end-1);
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
                   if tlist(j) == "(+)"
                       printTrace = "(Has GPR)";
                   elseif tlist(j) == "(X)"
                       printTrace = "(Reaction originated from master database)";
                   elseif tlist(j) == "(T)"
                       printTrace = "(Reaction is a transfer reaction)";
                   elseif tlist(j) == "(0)"
                       printTrace = "(Reaction originated from current model)";
                   elseif str2double(extractBetween(tlist(j),'(',')')) <= length(originData)
                       printTrace = "(" + originData(str2double(extractBetween(tlist(j),'(',')'))).origin + ")";
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
                       if length(groupList)>0
                           %groupSpacing = (blockR-blockL)/length(groupList)
                           groupSpacing = (blockR-blockL)/(length(groupList));
                           groupspace=blockL;
                           for j=1:length(groupList)
                               group = groupList(j);
                               gchr=char(groupList(j));
                               if endsWith(gchr,'x')
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
                checkMets_temp_removeLast = checkMets_temp(W,:)
                checkMets_temp(W,:) = checkMets_temp_removeLast(1:end-1);
                skippedMets = [skippedMets;initMet];

                dataStore_temp_edit{end+1, 1} = initMet;
                dataStore_temp_edit{end, 2} = 6;
                dataStore_temp_edit{end, 3} = trialStore;
                dataStore_temp_edit{end, 4} = pathLength;
                dataStore_temp_edit{end, 5} = randCnt_temp(W,1);
                dataStore_temp_edit{end, 6} = skipRandCnt_temp(W,1);
                dataStore_temp_edit{end, 7} = skipCnt_temp(W,1);
                dataStore_temp_edit{end, 8} = randSkipCnt_temp(W,1);
                dataStore_temp_edit{end, 9} = preRandRxnCnt;
                dataStore_temp_edit{end, 10} = randRxnCnt;
                dataStore_temp_edit{end, 11} = skipRxnCnt;
                dataStore_temp_edit{end, 12} = randSkipRxnCnt;
                dataStore_temp_edit{end, 13} = preRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 14} = randRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 15} = skipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 16} = randSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 17} = innerPreRandRxnCnt;
                dataStore_temp_edit{end, 18} = innerRandRxnCnt;
                dataStore_temp_edit{end, 19} = innerSkipRxnCnt;
                dataStore_temp_edit{end, 20} = innerRandSkipRxnCnt;
                dataStore_temp_edit{end, 21} = innerPreRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 22} = innerRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 23} = innerSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 24} = innerRandSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 25} = newMetsCnt;
                dataStore_temp_edit{end, 26} = randNewMetsCnt;
                dataStore_temp_edit{end, 27} = skipNewMetsCnt;
                dataStore_temp_edit{end, 28} = randSkipNewMetsCnt;
                dataStore_temp_edit{end, 29} = newMetsProdCnt;
                dataStore_temp_edit{end, 30} = skipNewMetsProdCnt;
                dataStore_temp_edit{end, 31} = innerNewMetsCnt;
                dataStore_temp_edit{end, 32} = innerRandNewMetsCnt;
                dataStore_temp_edit{end, 33} = innerSkipNewMetsCnt;
                dataStore_temp_edit{end, 34} = innerRandSkipNewMetsCnt;
                dataStore_temp_edit{end, 35} = cscore;
                dataStore_temp_edit{end, 36} = outerCScore;
                dataStore_temp_edit{end, 37} = outerCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 38} = randOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 39} = skipOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 40} = randSkipOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 41} = preRandCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 42} = randCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 43} = skipCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 44} = randSkipCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 45} = innerLoopCnt;
                dataStore_temp_edit{end, 46} = innerLoopCntTot_temp(W,1);
                dataStore_temp_edit{end, 47} = innerForLoopCntTot_temp(W,1);
                dataStore_temp_edit{end, 48} = outerLoopCnt_temp(W,1);
                dataStore_temp_edit{end, 49} = randSkipActive;
                dataStore_temp_edit{end, 50} = trialSuccess;

                fprintf(dataFile, '%s', initMet);
                fprintf(dataFile, '\t%d', dataStore_temp_edit{end, 2});
                fprintf(dataFile, '\t%d', trialStore);
                fprintf(dataFile, '\t%d', pathLength);
                fprintf(dataFile, '\t%d', randCnt_temp(W,1));
                fprintf(dataFile, '\t%d', skipRandCnt_temp(W,1));
                fprintf(dataFile, '\t%d', skipCnt_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipCnt_temp(W,1));
                fprintf(dataFile, '\t%d', preRandRxnCnt);
                fprintf(dataFile, '\t%d', randRxnCnt);
                fprintf(dataFile, '\t%d', skipRxnCnt);
                fprintf(dataFile, '\t%d', randSkipRxnCnt);
                fprintf(dataFile, '\t%d', preRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', randRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerPreRandRxnCnt);
                fprintf(dataFile, '\t%d', innerRandRxnCnt);
                fprintf(dataFile, '\t%d', innerSkipRxnCnt);
                fprintf(dataFile, '\t%d', innerRandSkipRxnCnt);
                fprintf(dataFile, '\t%d', innerPreRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerRandSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', newMetsCnt);
                fprintf(dataFile, '\t%d', randNewMetsCnt);
                fprintf(dataFile, '\t%d', skipNewMetsCnt);
                fprintf(dataFile, '\t%d', randSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', newMetsProdCnt);
                fprintf(dataFile, '\t%d', skipNewMetsProdCnt);
                fprintf(dataFile, '\t%d', innerNewMetsCnt);
                fprintf(dataFile, '\t%d', innerRandNewMetsCnt);
                fprintf(dataFile, '\t%d', innerSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', innerRandSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', cscore);
                fprintf(dataFile, '\t%d', outerCScore);
                fprintf(dataFile, '\t%d', outerCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', preRandCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerLoopCnt);
                fprintf(dataFile, '\t%d', innerLoopCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerForLoopCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', outerLoopCnt_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipActive);
                fprintf(dataFile, '\t%d\n', trialSuccess);

                %Resetting appropriate values
                randCnt_temp(W,1) = 0;
                skipRandCnt_temp(W,1) = 0;
                skipCnt_temp(W,1) = 0;
                randSkipCnt_temp(W,1) = 0;
                preRandRxnCnt = 0;
                randRxnCnt = 0;
                skipRxnCnt = 0;
                randskipRxnCnt = 0;
                preRandRxnCntTot_temp(W,1) = 0;
                randRxnCntTot_temp(W,1) = 0;
                skipRxnCntTot_temp(W,1) = 0; 
                randSkipRxnCntTot_temp(W,1) = 0;
                innerPreRandRxnCnt = 0;
                innerrandRxnCnt = 0; 
                innerskipRxnCnt = 0;
                innerRandSkipRxnCnt = 0;
                innerPreRandRxnCntTot_temp(W,1) = 0;
                innerRandRxnCntTot_temp(W,1) = 0;
                innerSkipRxnCntTot_temp(W,1) = 0;
                innerRandSkipRxnCntTot_temp(W,1) = 0;
                newMetsCnt = 0; 
                randnewMetsCnt = 0; 
                skipnewMetsCnt = 0;
                randSkipnewMetsCnt = 0;
                newMetsProdCnt = 0;
                skipnewMetsProdCnt = 0;
                innernewMetsCnt = 0;
                innerRandnewMetsCnt = 0;
                innerSkipnewMetsCnt = 0;
                innerRandSkipnewMetsCnt = 0;
                outerCScore = 0;
                outerCScoreTot_temp(W,1) = 0;
                randOuterCScoreTot_temp(W,1) = 0;
                skipOuterCScoreTot_temp(W,1) = 0;
                randSkipOuterCScoreTot_temp(W,1) = 0;
                preRandCScoreTot_temp(W,1) = 0;
                randCScoreTot_temp(W,1) = 0;
                skipCScoreTot_temp(W,1) = 0;
                randSkipCScoreTot_temp(W,1) = 0;
                innerLoopCnt = 0;
                innerLoopCntTot_temp(W,1) = 0;
                innerForLoopCntTot_temp(W,1) = 0;
                outerLoopCnt_temp(W,1) = 0;
           else
                trialSuccess = 3; 
                checkMets_temp_removeLast = checkMets_temp(W,:)
                checkMets_temp(W,:) = checkMets_temp_removeLast(1:end-1);
                fprintf("No solution found for " + initMet + " skipping over and continuing to next\n\n");
                failedMets = [failedMets;initMet];
                
                dataStore_temp_edit{end+1, 1} = initMet;
                dataStore_temp_edit{end, 2} = 7;
                dataStore_temp_edit{end, 3} = trialStore;
                dataStore_temp_edit{end, 4} = pathLength;
                dataStore_temp_edit{end, 5} = randCnt_temp(W,1);
                dataStore_temp_edit{end, 6} = skipRandCnt_temp(W,1);
                dataStore_temp_edit{end, 7} = skipCnt_temp(W,1);
                dataStore_temp_edit{end, 8} = randSkipCnt_temp(W,1);
                dataStore_temp_edit{end, 9} = preRandRxnCnt;
                dataStore_temp_edit{end, 10} = randRxnCnt;
                dataStore_temp_edit{end, 11} = skipRxnCnt;
                dataStore_temp_edit{end, 12} = randSkipRxnCnt;
                dataStore_temp_edit{end, 13} = preRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 14} = randRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 15} = skipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 16} = randSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 17} = innerPreRandRxnCnt;
                dataStore_temp_edit{end, 18} = innerRandRxnCnt;
                dataStore_temp_edit{end, 19} = innerSkipRxnCnt;
                dataStore_temp_edit{end, 20} = innerRandSkipRxnCnt;
                dataStore_temp_edit{end, 21} = innerPreRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 22} = innerRandRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 23} = innerSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 24} = innerRandSkipRxnCntTot_temp(W,1);
                dataStore_temp_edit{end, 25} = newMetsCnt;
                dataStore_temp_edit{end, 26} = randNewMetsCnt;
                dataStore_temp_edit{end, 27} = skipNewMetsCnt;
                dataStore_temp_edit{end, 28} = randSkipNewMetsCnt;
                dataStore_temp_edit{end, 29} = newMetsProdCnt;
                dataStore_temp_edit{end, 30} = skipNewMetsProdCnt;
                dataStore_temp_edit{end, 31} = innerNewMetsCnt;
                dataStore_temp_edit{end, 32} = innerRandNewMetsCnt;
                dataStore_temp_edit{end, 33} = innerSkipNewMetsCnt;
                dataStore_temp_edit{end, 34} = innerRandSkipNewMetsCnt;
                dataStore_temp_edit{end, 35} = cscore;
                dataStore_temp_edit{end, 36} = outerCScore;
                dataStore_temp_edit{end, 37} = outerCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 38} = randOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 39} = skipOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 40} = randSkipOuterCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 41} = preRandCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 42} = randCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 43} = skipCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 44} = randSkipCScoreTot_temp(W,1);
                dataStore_temp_edit{end, 45} = innerLoopCnt;
                dataStore_temp_edit{end, 46} = innerLoopCntTot_temp(W,1);
                dataStore_temp_edit{end, 47} = innerForLoopCntTot_temp(W,1);
                dataStore_temp_edit{end, 48} = outerLoopCnt_temp(W,1);
                dataStore_temp_edit{end, 49} = randSkipActive;
                dataStore_temp_edit{end, 50} = trialSuccess;

                fprintf(dataFile, '%s', initMet);
                fprintf(dataFile, '\t%d', dataStore_temp_edit{end, 2});
                fprintf(dataFile, '\t%d', trialStore);
                fprintf(dataFile, '\t%d', pathLength);
                fprintf(dataFile, '\t%d', randCnt_temp(W,1));
                fprintf(dataFile, '\t%d',skipRandCnt_temp(W,1));
                fprintf(dataFile, '\t%d', skipCnt_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipCnt_temp(W,1));
                fprintf(dataFile, '\t%d', preRandRxnCnt);
                fprintf(dataFile, '\t%d', randRxnCnt);
                fprintf(dataFile, '\t%d', skipRxnCnt);
                fprintf(dataFile, '\t%d', randSkipRxnCnt);
                fprintf(dataFile, '\t%d', preRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', randRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerPreRandRxnCnt);
                fprintf(dataFile, '\t%d', innerRandRxnCnt);
                fprintf(dataFile, '\t%d', innerSkipRxnCnt);
                fprintf(dataFile, '\t%d', innerRandSkipRxnCnt);
                fprintf(dataFile, '\t%d', innerPreRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerRandRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerRandSkipRxnCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', newMetsCnt);
                fprintf(dataFile, '\t%d', randNewMetsCnt);
                fprintf(dataFile, '\t%d', skipNewMetsCnt);
                fprintf(dataFile, '\t%d', randSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', newMetsProdCnt);
                fprintf(dataFile, '\t%d', skipNewMetsProdCnt);
                fprintf(dataFile, '\t%d', innerNewMetsCnt);
                fprintf(dataFile, '\t%d', innerRandNewMetsCnt);
                fprintf(dataFile, '\t%d', innerSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', innerRandSkipNewMetsCnt);
                fprintf(dataFile, '\t%d', cscore);
                fprintf(dataFile, '\t%d', outerCScore);
                fprintf(dataFile, '\t%d', outerCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipOuterCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', preRandCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', skipCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipCScoreTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerLoopCnt);
                fprintf(dataFile, '\t%d', innerLoopCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', innerForLoopCntTot_temp(W,1));
                fprintf(dataFile, '\t%d', outerLoopCnt_temp(W,1));
                fprintf(dataFile, '\t%d', randSkipActive);
                fprintf(dataFile, '\t%d\n', trialSuccess);

                randCnt_temp(W,1) = 0;
                skipRandCnt_temp(W,1) = 0;
                skipCnt_temp(W,1) = 0;
                randSkipCnt_temp(W,1) = 0;
                preRandRxnCnt= 0;
                randRxnCnt = 0;
                skipRxnCnt = 0;
                randskipRxnCnt = 0;
                preRandRxnCntTot_temp(W,1) = 0;
                randRxnCntTot_temp(W,1)= 0;
                skipRxnCntTot_temp(W,1) = 0; 
                randSkipRxnCntTot_temp(W,1) = 0;
                innerPreRandRxnCnt = 0;
                innerrandRxnCnt = 0; 
                innerskipRxnCnt = 0;
                innerRandSkipRxnCnt = 0;
                innerPreRandRxnCntTot_temp(W,1) = 0;
                innerRandRxnCntTot_temp(W,1) = 0;
                innerSkipRxnCntTot_temp(W,1) = 0;
                innerRandSkipRxnCntTot_temp(W,1) = 0;
                newMetsCnt = 0; 
                randnewMetsCnt = 0; 
                skipnewMetsCnt = 0;
                randSkipnewMetsCnt = 0;
                newMetsProdCnt = 0;
                skipnewMetsProdCnt = 0;
                innernewMetsCnt = 0;
                innerRandnewMetsCnt = 0;
                innerSkipnewMetsCnt = 0;
                innerRandSkipnewMetsCnt = 0;
                outerCScore = 0;
                outerCScoreTot_temp(W,1) = 0;
                randOuterCScoreTot_temp(W,1) = 0;
                skipOuterCScoreTot_temp(W,1) = 0;
                randSkipOuterCScoreTot_temp(W,1) = 0;
                preRandCScoreTot_temp(W,1) = 0;
                randCScoreTot_temp(W,1) = 0;
                skipCScoreTot_temp(W,1) = 0;
                randSkipCScoreTot_temp(W,1) = 0;
                innerLoopCnt = 0;
                innerLoopCntTot_temp(W,1) = 0;
                innerForLoopCntTot_temp(W,1) = 0;
                outerLoopCnt_temp(W,1) = 0;
           end
           exTrialFlag_temp(W,1)=0;
           exConfFlag_temp(W,1)=0;
           masterFlag_temp(W,1)=0;
           if ~largePathPerm
                maxTrials_temp(W,1) = stdTrials;
                largePathFlag_temp(W,1)=0;
                pathLength = stdLenMax;
           else
               maxTrials_temp(W,1) = largePathTrials;
            end
        end
        checkMets_temp_removeX = checkMets_temp(W,:)
        checkMets_temp(W,:) = checkMets_temp_removeX(~strcmp(checkMets_temp_removeX,"X"));
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

        fprintf(trialFile, '%s', initMet);
        fprintf(trialFile, '\t%s', fullRxnIDs);
        fprintf(trialFile, '\t%s', outerRxnIDs);
        fprintf(trialFile, '\t%s', innerRxnIDs);
        fprintf(trialFile, '\t%s', fullTracers);
        fprintf(trialFile, '\t%s', outerTracers);
        fprintf(trialFile, '\t%s', innerTracers);
        fprintf(trialFile, '\t%s', fullRxns);
        fprintf(trialFile, '\t%s', outerRxns);
        fprintf(trialFile, '\t%s', innerRxns);
        fprintf(trialFile, '\t%s', innerAllRxns);
        fprintf(trialFile, '\t%s', fullRevs);
        fprintf(trialFile, '\t%s', outerRevs);
        fprintf(trialFile, '\t%s', innerRevs);
        fprintf(trialFile, '\t%s', fullComps);
        fprintf(trialFile, '\t%s', outerComps);
        fprintf(trialFile, '\t%s', innerComps);
        fprintf(trialFile, '\t%s', fullTraceNums);
        fprintf(trialFile, '\t%s', outerTraceNums);
        fprintf(trialFile, '\t%s', innerTraceNums);
        fprintf(trialFile, '\t%d', ruleBreakLvl);
        fprintf(trialFile, '\t%d', skipAddActive);
        fprintf(trialFile, '\t%d', randPermActive);
        fprintf(trialFile, '\t%d', trialSuccess);
        fprintf(trialFile, '\t%d\n', failConfCheck);

%         trialInfo_temp_edit = trialInfo_temp{W,:,:}
%         trialInfo_temp_edit{end,2}  = fullRxnIDs;
%         trialInfo_temp_edit{end,3}  = outerRxnIDs;
%         trialInfo_temp_edit{end,4}  = innerRxnIDs;
%         trialInfo_temp_edit{end,5}  = fullTracers;
%         trialInfo_temp_edit{end,6}  = outerTracers;
%         trialInfo_temp_edit{end,7}  = innerTracers;
%         trialInfo_temp_edit{end,8}  = fullRxns;
%         trialInfo_temp_edit{end,9}  = outerRxns;
%         trialInfo_temp_edit{end,10} = innerRxns;
%         trialInfo_temp_edit{end,11} = innerAllRxns;
%         trialInfo_temp_edit{end,12} = fullRevs;
%         trialInfo_temp_edit{end,13} = outerRevs;
%         trialInfo_temp_edit{end,14} = innerRevs;
%         trialInfo_temp_edit{end,15} = fullComps;
%         trialInfo_temp_edit{end,16} = outerComps;
%         trialInfo_temp_edit{end,17} = innerComps;
%         trialInfo_temp_edit{end,18} = fullTraceNums;
%         trialInfo_temp_edit{end,19} = outerTraceNums;
%         trialInfo_temp_edit{end,20} = innerTraceNums;
%         trialInfo_temp_edit{end,21} = ruleBreakLvl;
%         trialInfo_temp_edit{end,22} = skipAddActive;
%         trialInfo_temp_edit{end,23} = randPermActive;
%         trialInfo_temp_edit{end,24} = trialSuccess;
%         trialInfo_temp_edit{end,25} = failConfCheck;
    %dataStore_temp{W} = dataStore_temp_edit
    end
%     trialInfo{W,:,:} = trialInfo_temp
    dataStore{W,:,:} = dataStore_temp_edit
    blockedMets_edit = blockedMets(W,:)
    blockedMets(W,:) = blockedMets_edit
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
fclose(trialFile);
fclose(dataFile);
warning('STOP');
tEnd = toc(tStart);

tracerResults.ID2add = ID2add;
tracerResults.rxn2add = rxn2add;
tracerResults.rev2add = rev2add;
tracerResults.boostID2add = boostID2add;
tracerResults.boostRxn2add = boostRxn2add;
tracerResults.boostRev2add = boostRev2add;
end

        



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    