function [addRxns] = tracerLoop2(met,traceNum,viableMets,excludeRxns,presentDB,superDB,orthoDB,masterDB,largePathFlag,masterFlag,netRxnsMax,skipAdd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

presRxns = string(presentDB.Abbreviation);

addRxns = {}; %list of reactions whose products contain the metabolite (met) in question
[eRxns] = tracerExistReactExtract(met,presentDB,traceNum); %pulls out reactions already in the model that produce the met in question
skipFlag = 0;
%check for reactions already in model that are producing the metabolite of
%interest, skip mets whose production should already be feasible
if ~isempty(eRxns) && length(eRxns(:,1))>1 && ~skipAdd
    for i = 1:length(eRxns(:,1))
        rt = eRxns{i,2}; %reactants
        pt = eRxns{i,3}; %products
        rx = eRxns{i,1}; %reaction id
        if ~ismember(rx,excludeRxns) && ~all(ismember(pt,viableMets)) %prevents reversible reactions causing the function to get stuck in an infinite loop
            addRxns = vertcat(addRxns, eRxns(i,:));
            %addRxns = eRxns(i,:);
        end
    end
end
if ~skipFlag    
    %Check for reactions in supermodel of which all reactants are viable
%     warning('Look here');
    if isempty(addRxns)
        [sRxns] = tracerSuperReactExtract(met,superDB,orthoDB,traceNum,presRxns); %pulls out reactions in the constructed "super" C4 plant model that produce the met in question
        sprt = []; %products
        if ~isempty(sRxns)
            for i = 1:length(sRxns(:,1))
                try
                    rt = sRxns{i,2}; %reactants
                catch
                end
                rx = sRxns{i,1}; %reaction ids
                sprt = [sprt,sRxns{i,2}];
                try
                    if all(ismember(rt,viableMets)) && ~ismember(rx,excludeRxns) %passes if all reactants of the reaction are viable and the reaction is not on the list of reactions to exclude
                        gene = sRxns{i,8};
                        [geneInfo] = checkGPRmatch(gene, orthoDB);  %Determining if rxn has a GPR match
                        if ~isempty(geneInfo{1})
                            sRxns{i, 6} = "(+)";
                        end
                        sRxns{i,8} = geneInfo{1};
                        sRxns{i,9} = geneInfo{2};
                        sRxns{i,10} = geneInfo{3};
                        addRxns = vertcat(addRxns,sRxns(i,:));
                    elseif largePathFlag && length(rt)-sum(ismember(rt,viableMets))<=netRxnsMax %alternatively will pass if larger pathways are allowed
                        gene = sRxns{i,8};
                        [geneInfo] = checkGPRmatch(gene, orthoDB);
                        if ~isempty(geneInfo{1})
                            sRxns{i, 6} = "(+)";
                        end
                        sRxns{i,8} = geneInfo{1};
                        sRxns{i,9} = geneInfo{2};
                        sRxns{i,10} = geneInfo{3};
                        addRxns = vertcat(addRxns,sRxns(i,:));
                    end
                catch
                end
            end
        end
    end
    if isempty(addRxns) && ~largePathFlag && ~isempty(sRxns) %if no qualifying reactions are found:
        target = string(mode(categorical((sprt)))); %sets target metabolite to the product showing up the most times in super model reaction reactants
        idxList = find(contains(string(sRxns(:,4)),target)); %pulls out reactions producing the target
        rxcov = [];
        for i = 1:length(idxList)
            idx = idxList(i);
            rt = string(sRxns{idx,2});
            sum(ismember(rt,viableMets)/length(rt));
            rxcov = [rxcov;(sum(ismember(rt,viableMets))/length(rt))]; %calculates how many viable metabolites as products for each rxn
        end
        fidx = find(rxcov==max(rxcov)); %finds index of reaction with most viable metabolites as products
        for i = 1:length(fidx)
            rx = sRxns{fidx(i),1}; %finds associated reaction ids
            if  ~ismember(rx,excludeRxns) %checks that the found reaction is not in the excluded reaction list
                gene = sRxns{fidx(i),8};
                [geneInfo] = checkGPRmatch(gene, orthoDB);
                if ~isempty(geneInfo{1})
                    sRxns{fidx(i), 6} = "(+)";
                end
                sRxns{fidx(i),8} = geneInfo{1};
                sRxns{fidx(i),9} = geneInfo{2};
                sRxns{fidx(i),10} = geneInfo{3};
                addRxns = vertcat(addRxns, sRxns(fidx(i), :));
                %addRxns = sRxns(fidx(i),:); %adds reaction as a potential part of solution
            end
        end
        if isempty(addRxns) %if no reaction was found:
            for i = 1:length(idxList)
               idx = idxList(i);
                rx = sRxns{idx,1};
                nrt = sRxns{idx,2};
                if ~ismember(rx,excludeRxns) && length(rt)-sum(ismember(rt,viableMets))<=1 %adds a reaction that could potentially solve but is not a great candidate becuase multiple products are missing
                    gene = sRxns{idx,8};
                    [geneInfo] = checkGPRmatch(gene, orthoDB);
                    if ~isempty(geneInfo{1})
                        sRxns{idx, 6} = "(+)";
                    end
                    sRxns{idx,8} = geneInfo{1};
                    sRxns{idx,9} = geneInfo{2};
                    sRxns{idx,10} = geneInfo{3};
                    addRxns = vertcat(addRxns, sRxns(idx, :));
                    %addRxns = sRxns(idx,:);
                end
            end
        end
    end
    if isempty(addRxns)
        %Determine if transfer reaction is needed, transport from cytosol if
        %possible, if found in one other organelle other than cytosol,
        %transport to cytosol first, if found in more than one organelle and
        %NOT cytosol,pass and find if identical reactions producing metabolite
        %are occuring in each organelle, if yes, add reaction to organelle and
        %cytosol
        if contains(met,"_mr") 
            metspl = strsplit(met,"_");
            mettag = "_"+metspl(end);
            mcomp = extractBetween(mettag,"[","]"); %stores the compartment code ex. "_[c]"
            nakedMet = metspl(1:end-1); %naked met is metabolite id with the compartment code stripped off ex. "C00001"
            vmets = viableMets(find(contains(viableMets,nakedMet))); %finds viable metabolites in all compartments matching the one in question
            vstr = string(vmets);
            vspl = arrayfun(@(s) strsplit(s,"_"),vstr,'UniformOutput',false); %seperates metabolite ids from the compartment code
            vtag = [];
            vcomp = [];
            for i = 1:length(vspl)
                vt = vspl{i}(2);
                vtag = [vtag;vt];
                vcomp = [vcomp;extractBetween(vt,"[","]")]; %removes formating from compartment code ex. "_[c]" --> "c"
            end
        else
            metspl = strsplit(met,{'[',']'});
            mettag = "["+metspl{2}+"]";
            mcomp = metspl{2}; %stores the compartment code ex. "_[c]"
            nakedMet = metspl{1}; %naked met is metabolite id with the compartment code stripped off ex. "C00001"
            vmets = viableMets(find(contains(viableMets,nakedMet))); %finds viable metabolites in all compartments matching the one in question
            vstr = string(vmets);
            vspl = arrayfun(@(s) strsplit(s,{'[',']'}),vstr,'UniformOutput',false); %seperates metabolite ids from the compartment code
            vtag = [];
            vcomp = [];
            for i = 1:length(vspl)
                vt = "[" + vspl{i}(2) + "]";
                vtag = [vtag;vt];
                vcomp = [vcomp;vspl{i}(2)]; %removes formating from compartment code ex. "_[c]" --> "c"
            end
        end
        
        if ~isempty(vcomp)
            vcnt = sum(~contains(vcomp,"c"));
            vcompExl = vcomp(~contains(vcomp,"c"));
            vtagExl = vtag(~contains(vtag,"[c]"));
            if any(contains(vcomp,"c")) && mcomp ~= "c" 
                cytotag = strrep(mettag,"["+mcomp+"]","[c]");
                transID = "Transport_c"+mcomp+"_"+met;
                if ~ismember(transID,excludeRxns)
                    addRxns = {transID,nakedMet+cytotag,met,nakedMet+cytotag+" <=> "+met,1,"(T)",traceNum};
                end
            elseif mcomp == "c" && vcnt == 1
                transtag = strrep(mettag,"[c]","["+vcompExl+"]");
                transID = "Transport_c"+vcompExl+"_"+nakedMet+transtag;
                if ~ismember(transID,excludeRxns)
                    addRxns = {transID,nakedMet+transtag,met,met+" <=> "+nakedMet+transtag,1,"(T)",traceNum};
                end
            elseif mcomp ~= "c" && vcnt == 1
                cytotag = strrep(mettag,"["+mcomp+"]","[c]");
                transtag = strrep(mettag,"["+mcomp+"]","["+vcompExl+"]");
                transID1 = "Transport_c"+vcompExl+"_"+nakedMet+transtag;
                outRxn1 = {transID1,nakedMet+transtag,nakedMet+cytotag,nakedMet+cytotag+" <=> "+nakedMet+transtag,1,"(T)",traceNum};
                transID2 = "Transport_c"+mcomp+"_"+met;
                outRxn2 = {transID2,nakedMet+cytotag,met,nakedMet+cytotag+" <=> "+met,1,"(T)",traceNum};
                if ~ismember(transID2,excludeRxns) && ~ismember(transID1,excludeRxns)
                    addRxns = vertcat(outRxn1,outRxn2);
                end
            end
        end
    end
    if ~iscell(addRxns)
        addRxns = {};
    % search master reaction database for a reaction that fits. This operates similarly the super reaction processing above   
    elseif masterFlag
        [mRxns] = tracerMasterReactExtract(met,masterDB,traceNum, presRxns);
        if ~isempty(mRxns)
            parfor i = 1:length(mRxns(:,1)),5
                rt = mRxns{i,2};
                rx = mRxns{i,1};
                if all(ismember(rt,viableMets)) && ~any(strcmp(rx,excludeRxns))
                    addRxns = vertcat(addRxns,mRxns(i,:));
                elseif largePathFlag && length(rt)-sum(ismember(rt,viableMets))<=netRxnsMax
                    addRxns = vertcat(addRxns,mRxns(i,:));
                end
            end
        end
    end
else 
    addRxns = "SKIP"; %assigns dummy "SKIP" string signifying that it metabolite is not being produced for reasons not relating to a missing reaction
end

function[geneInfo] = checkGPRmatch(gene, orthoDB);
zmaOr = orthoDB.zma;
sbiOr = orthoDB.sbi;
araOr = orthoDB.ara;
sitaOr = orthoDB.sita;
bdiOr = orthoDB.bdi;
targetSpeciesOr = zmaOr;
alphaNum = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z' '0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};
num = {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};
alpha = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
geneBySpecies = extractBetween(gene, "[", "]");
            

    %Extracts each gene option for each species (each cell in a row
    %is a gene or set of genes that work together

    %Determining if each gene has a known zma gene
    %
    listGPR = '';                       %stores all GPRs with orthologs in target species
    listOrtho = '';                     %stores all orthologs in target species that map to GPR
    listSpecies = '';                   %Stores all species origins for each gene with an ortholog match

    %Loops through each gene grouping
    for speciesIdx = 1:length(geneBySpecies) 
        %Extracts each gene or set of genes that must work together by parentheses "( )"(linked by &'s)
        geneOptions = extractBetween(geneBySpecies(speciesIdx), "(", ")"); 
        %Loops through each gene in grouping
        for geneIdx = 1:length(geneOptions) 
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
                   % - zma
                    if any(contains(zmaOr, genespl)) 
                        idxOr = find(contains(zmaOr, genespl));
                        orthoSpecies = 'zma';
                    % - sbi
                    elseif any(contains(sbiOr, genespl))
                        idxOr = find(contains(sbiOr, genespl));
                        orthoSpecies = 'sbi';
                    % - sita
                    elseif any(contains(sitaOr, genespl))
                        idxOr = find(contains(sitaOr, genespl));
                        orthoSpecies = 'sita';
                    % - ara
                    elseif any(contains(araOr, genespl))
                        idxOr = find(contains(araOr, genespl));
                        orthoSpecies = 'ara';
                    % - bdi
                    elseif any(contains(bdiOr, genespl))
                        idxOr = find(contains(bdiOr, genespl));
                        orthoSpecies = 'bdi';
                    end 

                    %Initialize cell to store the gene IDs of the
                    %target species's genes that they match with
                    targetOrthoStore = {};

                    %If one or more matches were returned,
                    %determine if a target species gene exists for the gene
                    if sum(idxOr) ~= 0      
                        for k = 1:length(idxOr)
                            %extracting next gene in target species
                            %that matches
                            targetOrtho = char(targetSpeciesOr(idxOr(k)));

                            %Determines if there is a known,
                            %corresponding gene in target species

                            %If the target species gene exists and
                            %has not already been stored
                            if ~isempty(targetOrthoStore) & ~any(contains(targetOrthoStore, targetOrtho)) & contains(targetOrtho, alphaNum, 'IgnoreCase', true) %This might need to be changed to ~isnan
                                targetOrthoStore{end+1} = targetOrtho;  %Stores the ortholog in target species so repeats are not added to the output
                                rank="+";                               %Assigning GPR match rank

                            elseif isempty(targetOrthoStore) & contains(targetOrtho, alphaNum,'IgnoreCase', true)
                                targetOrthoStore{1} = targetOrtho;
                                rank = "+";                             %Assigning GPR match rank

                            end
                        end
                    end
                    %If a matching gene is found in target species,
                    %store it
                    if ~isempty(targetOrthoStore)
                        %Storing target species's orthologs for
                        %later use (if GPR match is found)
                        targetOrthoForm = ''; %array stores properly formatted possible genes for use 
                        for k = 1:length(targetOrthoStore)
                            targetOrthoForm = [targetOrthoForm, targetOrthoStore{k}, ' or '];
                        end
                        targetOrthoPass = targetOrthoForm(1:end-4); %Removing final 'or' operator
                        targetOrthoOut = ['[', targetOrthoPass, ']']; %Enclosing list in square brackets
                    end
                %end
            elseif contains(genespl, alpha, 'IgnoreCase', true)     %If this reaction requires multiple genes (there is an & in the GPR)
                geneComplex = strsplit(genespl, ' & ');             %Splits GPR into each gene--produces a cell array output
                %Determines if the gene exists in orthoDB, 
                %Determining what species the gene belongs to
                   %and finding indexes of matches
                   %orthoSpecies variable allows function to keep
                   %track of where ortholog exists
                % -zma
                if all(ismember(geneComplex, zmaOr))
                    andOrtho = zmaOr;
                    orthoSpecies = 'zma';
                % -sbi
                elseif all(ismember(geneComplex, sbiOr))
                    andOrtho = sbiOr;
                    orthoSpecies = 'sbi';
                % -sita
                elseif all(ismember(geneComplex, sitaOr))
                    andOrtho = sitaOr;
                    orthoSpecies = 'sita';
                % -ara
                elseif all(ismember(geneComplex, araOr))
                    andOrtho = araOr;
                    orthoSpecies = 'ara';
                % -bdi
                elseif all(ismember(geneComplex, bdiOr))
                    andOrtho = bdiOr;
                    orthoSpecies = 'bdi';
                else
                    andOrtho = '';
                end
                %Loop determines if each gene has an ortholog in
                %target species and extracts those orthologs
                matchCount = 0;                                     %Stores how many genes also have an ortholog in the target species, after going through each gene, if this number is equal to the number of genes link by the '&' operator, then there is a GPR match
                targetOrthoMatch = '';
                for andIdx = 1:length(geneComplex)
                    andGene = geneComplex{andIdx};                  %Extracts each individual gene needed
                    idxOr = find(contains(andOrtho, andGene));      %Finds locations of matching gene in the appropriate ortholog gene list

                    %If there is 1+ ortholog in the target species
                    if contains(targetSpeciesOr(idxOr), alphaNum, 'IgnoreCase', true)
                        targetOrthoStore = {}; %Cell array stores all target species orthologs for given gene
                        %Loop to find target species orthologs for
                        %current gene
                        matchCount = matchCount + 1;                %Stores how many of the required ortholog genes have a match in the target species
                        for retIdx = 1:length(idxOr)
                            targetOrtho = char(targetSpeciesOr(idxOr(retIdx))); %Retrieves next gene in target species that is an ortholog to the current gene
                            %If there is a gene and it has not
                            %already been stored:
                            if ~isempty(targetOrthoStore) & ~any(contains(targetOrthoStore, targetOrtho)) & contains(targetOrtho, alphaNum, 'IgnoreCase', true)
                                targetOrthoStore{end+1} = targetOrtho;                          %Adding target species ortholog to list of orthologs so repeats are not added
                            %If this is the first gene stored
                            elseif isempty(targetOrthoStore) & contains(targetOrtho, alphaNum, 'IgnoreCase', true)
                                targetOrthoStore{1} = targetOrtho; 
                            end
                        end
                        %Storing target species's orthologs for later use (if GPR match is found)
                        targetOrthoForm = ''; %array stores properly formatted possible genes for use 
                        for k = 1:length(targetOrthoStore)
                            targetOrthoForm = [targetOrthoForm, targetOrthoStore{k}, ' or ']; %List of genes in target species that are known orthologs for current gene, now separated by 'or'
                        end
                        targetOrthoForm = ['(', targetOrthoForm(1:end-4),')']; %List of genes in target species that are known orthologs for current gene, now enclosed by parentheses
                        targetOrthoMatch = [targetOrthoMatch, targetOrthoForm, ' and ']; %List of all possible orthologs in target species that satisfy the GPR
                    end
                end
                %If each gene linked by '&' has an ortholog in the
                %target species
                if matchCount == length(geneComplex)
                    rank = "+"; %Assigning GPR rank to reaction
                    %Storing target species orthologs for output
                    targetOrthoPass = targetOrthoMatch(1:end-5); %Removing final 'and' from string
                    targetOrthoOut = ['[', targetOrthoPass, ']']; %Enclosing list in square brackets
                end
            end
            %Storing GPR, orthologs, and orthologs's species of
            %origin
           if ~isempty(targetOrthoOut)
                listGPR = [listGPR,'[',genespl,']'];
                listOrtho = [listOrtho, targetOrthoOut];
                listSpecies = [listSpecies, '[',orthoSpecies,']'];
           elseif speciesIdx == length(geneBySpecies) && isempty(listOrtho)
               listGPR = [];
               listOrtho = [];
               listSpecies = [];
           end
        end
    end
geneInfo{1} = listGPR;
geneInfo{2} = listOrtho;
geneInfo{3} = listSpecies;
end
end
