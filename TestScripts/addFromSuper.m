function [addRxns] = addFromSuper(met,traceNum,viableMets,excludeRxns,superDB,orthoDB, presentDB)

addRxns = {}; %list of reactions whose products contain the metabolite (met) in question
presRxns = string(presentDB.Abbreviation);



%Check for reactions in supermodel of which all reactants are viable
[sRxns] = tracerSuperReactExtract(met,superDB,orthoDB,traceNum,presRxns); %pulls out reactions in the constructed "super" C4 plant model that produce the met in question
sprt = []; %products
if ~isempty(sRxns)
    for i = 1:length(sRxns(:,1))
        try
            rt = sRxns{i,2}; %reactants
            pt = sRxns{i,3}; %products
            rev = sRxns{i,5}; %reversibility
        catch
        end
        rx = sRxns{i,1}; %reaction ids
        sprt = [sprt,sRxns{i,2}];
        if all(ismember(rt,viableMets)) && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) %passes if all reactants of the reaction are viable and the reaction is not on the list of reactions to exclude
            gene = sRxns{i,8};
            nG = 0;
            [nG] = checkGPRmatch(gene, orthoDB);  %Determining if rxn has a GPR match
            rank = str2double(extractBetween(sRxns{i,6},"(",")"));
            cscore = calcConfidenceScore(nG, rank, orthoDB, superDB);
            sRxns{i, 9} = cscore;  
            sRxns{i, 10} = 0;
            sRxns{i,17} = "None";
            addRxns = vertcat(addRxns,sRxns(i,:));
            
        elseif all(ismember(pt,viableMets)) && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) && rev == 1
            gene = sRxns{i,8};
            nG = 0;
            [nG] = checkGPRmatch(gene, orthoDB);  %Determining if rxn has a GPR match
            rank = str2double(extractBetween(sRxns{i,6},"(",")"));
            cscore = calcConfidenceScore(nG, rank, orthoDB, superDB);
            sRxns{i, 9} = cscore;  
            sRxns{i, 10} = 0;
            sRxns{i,17} = "None";
            addRxns = vertcat(addRxns,sRxns(i,:));
            
        elseif length(rt) - sum(ismember(rt,viableMets)) == 1 && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) %reactions where only ONE reactant is missing - attempts to solve will require multiple reactions
            gene = sRxns{i,8};
            nG = 0;
            [nG] = checkGPRmatch(gene, orthoDB);  %Determining if rxn has a GPR match
            rank = str2double(extractBetween(sRxns{i,6},"(",")"));
            cscore = calcConfidenceScore(nG, rank, orthoDB, superDB);
            sRxns{i, 9} = cscore;  
            sRxns{i, 10} = 1;
            missingMet = rt(~ismember(rt,viableMets));
            sRxns{i,17} = missingMet;
            addRxns = vertcat(addRxns,sRxns(i,:));
        elseif length(rt) - sum(ismember(pt,viableMets)) == 1 && ~any(ismember(pt,excludeRxns)) && ~any(ismember(rx,presRxns)) && rev == 1 %reactions where only ONE reactant is missing - attempts to solve will require multiple reactions
            gene = sRxns{i,8};
            nG = 0;
            [nG] = checkGPRmatch(gene, orthoDB);  %Determining if rxn has a GPR match
            rank = str2double(extractBetween(sRxns{i,6},"(",")"));
            cscore = calcConfidenceScore(nG, rank, orthoDB, superDB);
            sRxns{i, 9} = cscore;  
            sRxns{i, 10} = 1;
            missingMet = pt(~ismember(pt,viableMets));
            sRxns{i,17} = missingMet;
            addRxns = vertcat(addRxns,sRxns(i,:));            
        end
       
    end
end

try
    addRxns = sortrows(addRxns,[10,9],['descend','ascend']);
catch
end
end

function [cscore] = calcConfidenceScore(nG, rank,orthoDB,superDB)
speciesNames = orthoDB.Properties.VariableNames;
numOrthoSpecies = length(speciesNames);
maxRank = max(superDB.Rank);
pWeight = 0.5;
gWeight = 0.5;

if string(nG) == "MATCH"
    cscore = 1.0;
else
    gFactor = nG/numOrthoSpecies;
    pFactor = rank/maxRank;
    cscore = pWeight*pFactor + gWeight*gFactor;
end
end

function[nG] = checkGPRmatch(gene, orthoDB)
speciesNames = orthoDB.Properties.VariableNames;
orthoArray = orthoDB.Variables;
nG = 0;
nG_list = [];
% zmaOr = orthoDB.zma;
% sbiOr = orthoDB.sbi;
% araOr = orthoDB.ara;
% sitaOr = orthoDB.sita;
% bdiOr = orthoDB.bdi;
%targetSpeciesOr = zmaOr;
alphaNum = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z' '0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};
num = {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};
alpha = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
geneBySpecies = extractBetween(gene, "[", "]");
            
    [~, colsDB] = size(orthoDB);
    %Extracts each gene option for each species (each cell in a row
    %is a gene or set of genes that work together

    %Determining if each gene has a known zma gene
    %
    listGPR = '';                       %stores all GPRs with orthologs in target species
    listOrtho = '';                     %stores all orthologs in target species that map to GPR
    listSpecies = '';                   %Stores all species origins for each gene with an ortholog match
    nG_count =  zeros(1,colsDB);      %Creates an array to determine how many species have an ortholog to the gene in model species
    nG_counter = 0;
    %Loops through each gene grouping
    nG_list = [];
    for speciesIdx = 1:length(geneBySpecies) 
        %Extracts each gene or set of genes that must work together by parentheses "( )"(linked by &'s)
        geneOptions = extractBetween(geneBySpecies(speciesIdx), "(", ")"); 
        %Loops through each gene in grouping
        for geneIdx = 1:length(geneOptions)
            nG = 0;
            targetOrthoOut = '';                                    %Stores all orthologs for GPR
            targetOrthoResults = '';                                %Stores all orthologs for all GPRS (returned to main function list results)
            genespl = geneOptions{geneIdx};                         %Current gne of interest
            idxOr = [];                                             %Stores locations of all orthologs that are also located in orthoDB


        %                     This section determines if the genes have a maize
        %                     ortholog -> allowing for a GPR match to be determined:
        %                     --If the gene does not require other genes and is just a single gene
            if ~contains(genespl, '&') && contains(genespl, alpha, 'IgnoreCase', true)     
                %Determining what species the gene belongs to
                %and finding indexes of matches
                %orthoSpecies variable allows function to keep
                %track of where ortholog exists
                foundFlag = 0;
                colIdx = 1;
                numSpecies = length(speciesNames);
                emptyCells = {};
                
                for species_idx = 1:numSpecies
                    idxOr = find(strcmp(genespl, orthoArray(:,species_idx)));
                    if ~isempty(idxOr)
                        numOrth = length(idxOr);
                        %orthMatchTotals = zeros(numOrth,1);
                        orth2check = zeros(numOrth,1);
                        for or_idx = 1:numOrth
                            check_idx = idxOr(or_idx);
                            if ~isempty(orthoArray(check_idx, 1))
                                nG_list = "MATCH";
                            end
                            numMatch = 0;
                            %disp(orthoArray(check_idx,:))
                            emptyCells = cellfun('isempty',orthoArray(check_idx,:));
                            numMatch = length(emptyCells) - sum(emptyCells);
                            %orthMatchTotals(or_idx,1) = numMatch;
                            nG_list = [nG_list;numMatch];
                            nG_counter = nG_counter+1;
                            %disp("NUMMATCH")
                            %disp(numMatch)
                        end
                    else
                        nG_list = [nG_list;0];
                        nG_counter = nG_counter+1;
                    end
                end          
            end   
        end        
    end
if any(strcmp(nG_list,"MATCH"))
    nG = "MATCH";
elseif isempty(nG_list)
    nG = 0;
else
    nG = max(nG_list);
end
end

function [outRxns] = tracerSuperReactExtract(met,superDB,orthoDB,traceNum,presRxns)
keggFlag = 0;
supRxns = string(superDB.Reaction);
supRanks = superDB.Rank;
supSid = string(superDB.Abbreviation);
supKid = string(superDB.KeggID);
supGene = string(superDB.GPR);
try
    supLB = superDB.Lowerbound;
    supUB = superDB.Upperbound;
catch
    supLB = superDB.LowerBound;
    supUB = superDB.UpperBound;
end
supDescrip = string(superDB.Description);
supDef = string(superDB.Definition);
supDelG = superDB.DeltaG;
supDelGerror = superDB.DeltaGError;
supEC = string(superDB.ECnumber);
supPath = string(superDB.Pathway);

%ortholog lists
zmaOr = orthoDB.zma;
sbiOr = orthoDB.sbi;
araOr = orthoDB.ara;
sitaOr = orthoDB.sita;
bdiOr = orthoDB.bdi;
targetSpeciesOr = zmaOr;
alphaNum = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z' '0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};
num = {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};
alpha = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};

outRxns = {};
%strip tag from metabolite
if contains(met, "_mr")
    metspl = strsplit(met,"_");
    nakedMet = metspl(1);
    tag = "_"+metspl(2);
else
    metspl = strsplit(met,{'[',']'});
    nakedMet = metspl(1);
    tag = "[" + metspl(2) + "]";
end
    
%determine rxns producing met of interest
idxE = find(contains(supRxns,nakedMet));
if ~isempty(idxE)
    cnt = 0;
    for i = 1:length(idxE)
        idx1 = idxE(i);
        clear rxnRecon;
        pdt = [];
        rct = [];
        rank = supRanks(idx1);
        rxn = supRxns(idx1);
        sid = supSid(idx1);
        descrip = supDescrip(idx1);
        def = supDef(idx1);
        delG = supDelG(idx1);
        delGerror = supDelGerror(idx1);
        ec = supEC(idx1);
        path = supPath(idx1);
        if ~any(contains(presRxns, sid))% && rank ~= 1
            kid = supKid(idx1);
            ub = supUB(idx1);
            lb = supLB(idx1);
            gene = supGene(idx1);
            if lb < -1 && ub > 1
                rev = 1;
            else
                rev = 0;
            end
            if contains(rxn," <=> ")
                spl = strip(strsplit(rxn," <=> "));
                Lstr = spl(1);
                Rstr = spl(2);
                dir = " <=> ";
            elseif contains(rxn," => ") || contains(rxn," -> ")
                spl = strip(strsplit(rxn,"=>"));
                Lstr = spl(1);
                Rstr = spl(2);
                dir = " -> ";
            elseif contains(rxn," <= ") || contains(rxn," <- ")
                %Not the case for utilized database, consider ammending if time
            else
                warning("Master Database has unacceptable format for reaction '" + sid + "'")
                continue
            end
            flipFlag = 0;
            if lb < 0 && contains(Lstr,nakedMet)
                flipFlag = 1;
                prodList = strip(strsplit(Lstr," "));
                Lrecon = strings(size(prodList));
                for idx2 = 1:length(prodList)
                    prodstr = prodList(idx2);
                    if isnan(str2double(prodstr)) && ~startsWith(prodstr,"+") && ~startsWith(prodstr,"=") && ~startsWith(prodstr,"<") && ~startsWith(prodstr,"-")
                        if keggFlag
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            pdt = [pdt,(prodstr+tag)];
                            Lrecon(idx2) = prodstr+tag;
                        end

                    else
                        Lrecon(idx2) = prodstr;
                    end
                end
                reactList = strip(strsplit(Rstr," "));
                Rrecon = strings(size(reactList));
                for idx2 = 1:length(reactList)
                    reactstr = reactList(idx2);
                    if isnan(str2double(reactstr)) && ~startsWith(reactstr,"+") && ~startsWith(reactstr,"=") && ~startsWith(reactstr,"<") && ~startsWith(reactstr,"-")
                        if keggFlag
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            rct = [rct,(reactstr+tag)];
                            Rrecon(idx2) = reactstr+tag;
                        end

                    else
                        Rrecon(idx2) = reactstr;
                    end
                end
            elseif ub > 0 && contains(Rstr,nakedMet)
                prodList = strip(strsplit(Rstr," "));
                Rrecon = strings(size(prodList));
                for idx2 = 1:length(prodList)
                    prodstr = prodList(idx2);
                    if isnan(str2double(prodstr)) && ~startsWith(prodstr,"+") && ~startsWith(prodstr,"=") && ~startsWith(prodstr,"<") && ~startsWith(prodstr,"-")
                        if keggFlag
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            pdt = [pdt,(prodstr+tag)];
                            Rrecon(idx2) = prodstr+tag;
                        end

                    else
                        Rrecon(idx2) = prodstr;
                    end
                end  
                reactList = strip(strsplit(Lstr," "));
                Lrecon = strings(size(reactList));
                for idx2 = 1:length(reactList)
                    reactstr = reactList(idx2);
                    if isnan(str2double(reactstr)) && ~startsWith(reactstr,"+") && ~startsWith(reactstr,"=") && ~startsWith(reactstr,"<") && ~startsWith(reactstr,"-")
                        if keggFlag
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            rct = [rct,(reactstr+tag)];
                            Lrecon(idx2) = reactstr+tag;
                        end

                    else
                        Lrecon(idx2) = reactstr;
                    end
                end
            else
                Lrecon = [];
                Rrecon = [];
            end                 
            if ~isempty(Lrecon) && ~isempty(Rrecon)
                if flipFlag
                    rxnRecon = strjoin(Rrecon," ") + dir + strjoin(Lrecon," ");
                else
                    rxnRecon = strjoin(Lrecon," ") + dir + strjoin(Rrecon," ");
                end
            elseif ~isempty(Lrecon) && isempty(Rrecon)
                warning("Reaction, " + rxn + " has no right-hand metabolites")
                rxnRecon = strip(strjoin(Lrecon," ") + dir);
            elseif isempty(Lrecon) && ~isempty(Rrecon)
                warning("Reaction, " + rxn + " has no left-hand metabolites")
                rxnRecon = strip(dir + strjoin(Rrecon," "));
            end
            %%
            if exist('rxnRecon','var')

                %Format output information for each rxn
                cnt = cnt+1;
                outRxns{cnt,1} = sid + tag;
                outRxns{cnt,2} = rct;
                outRxns{cnt,3} = pdt;
                outRxns{cnt,4} = strip(rxnRecon);
                outRxns{cnt,5} = rev;
                outRxns{cnt,6} = "("+rank+")";
                outRxns{cnt,7} = traceNum;
                outRxns{cnt,8} = gene;
                outRxns{cnt,11} = descrip;
                outRxns{cnt,12} = def;
                outRxns{cnt,13} = delG;
                outRxns{cnt,14} = delGerror;
                outRxns{cnt,15} = ec;
                outRxns{cnt,16} = path;
                
            end
        end
    end
end
end


