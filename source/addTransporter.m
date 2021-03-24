function [addRxns] = addTransporter(met,traceNum,viableMets,excludeRxns,transDB, presentDB)

addRxns = {}; %list of reactions whose products contain the metabolite (met) in question
presRxns = string(presentDB.Abbreviation);
compartments = string(["[c]", "[m]", "[x]", "[d]", "[r]", "[e]"]);


%Check for reactions in supermodel of which all reactants are viable
[tRxns] = transportReactExtract(met,transDB,traceNum,presRxns, viableMets); %pulls out reactions in the constructed "super" C4 plant model that produce the met in question
sprt = []; %products
if ~isempty(tRxns)
    for i = 1:length(tRxns(:,1))
        try
            rt = tRxns{i,2}; %reactants
            pt = tRxns{i,3}; %products
            rev = tRxns{i,5}; %reversibility
        catch
        end
        rx = tRxns{i,1}; %reaction ids
        sprt = [sprt,tRxns{i,2}];
        if all(ismember(rt,viableMets)) && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) %passes if all reactants of the reaction are viable and the reaction is not on the list of reactions to exclude
            rank = "(T)";
            cscore = 0;
            tRxns{i, 9} = cscore;  
            tRxns{i, 10} = 0;
            tRxns{i,17} = "None";
            addRxns = vertcat(addRxns,tRxns(i,:));
            
        elseif all(ismember(pt,viableMets)) && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) && rev == 1
            rank = "(T)";
            cscore = 0;
            tRxns{i, 9} = cscore;  
            tRxns{i, 10} = 0;
            tRxns{i,17} = "None";
            addRxns = vertcat(addRxns,tRxns(i,:));
            
        elseif length(rt) - sum(ismember(rt,viableMets)) == 1 && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) %reactions where only ONE reactant is missing - attempts to solve will require multiple reactions
            rank = "(T)";
            cscore = 0;
            tRxns{i, 9} = cscore;  
            tRxns{i, 10} = 1;
            missingMet = rt(~ismember(rt,viableMets));
            tRxns{i,17} = missingMet;
            addRxns = vertcat(addRxns,tRxns(i,:));
        elseif length(rt) - sum(ismember(pt,viableMets)) == 1 && ~any(ismember(pt,excludeRxns)) && ~any(ismember(rx,presRxns)) && rev == 1 %reactions where only ONE reactant is missing - attempts to solve will require multiple reactions
            rank = "(T)";
            cscore = 0;
            tRxns{i, 9} = cscore;  
            tRxns{i, 10} = 1;
            missingMet = pt(~ismember(pt,viableMets));
            tRxns{i,17} = missingMet;
            addRxns = vertcat(addRxns,tRxns(i,:));            
        end
       
    end
end

try
    addRxns = sortrows(addRxns,[10,9],['descend','ascend']);
catch
end
end


function [outRxns] = transportReactExtract(met,transDB,traceNum,presRxns,viableMets)
keggFlag = 0;
transRxns = string(transDB.Reaction);
transSid = string(transDB.Abbreviation);
try
    transLB = transDB.Lowerbound;
    transUB = transDB.Upperbound;
catch
    transLB = transDB.LowerBound;
    transUB = transDB.UpperBound;
end
transDescrip = string(transDB.Description);
transDef = string(transDB.Definition);
transDelG = transDB.DeltaG;
transDelGerror = transDB.DeltaGError;
transEC = string(transDB.ECnumber);

compartments = string(["[c]", "[m]", "[x]", "[d]", "[r]", "[e]"]);
outRxns = {};
%strip tag from metabolite
metspl = strsplit(met,"[");
nakedMet = metspl(1);
tag = "[" + extractBetween(met,"[","]") + "]";
%determine rxns producing met of interest
idxE = find(contains(transRxns,nakedMet));
if ~isempty(idxE)
    cnt = 0;
    for i = 1:length(idxE)
        for c = 1:length(compartments)
            comptag = compartments(c);
            compMet = nakedMet + comptag;
            if ~contains(viableMets,compMet)
                continue
            end
            idx1 = idxE(i);
            clear rxnRecon;
            pdt = [];
            rct = [];
            rxn = transRxns(idx1);
            sid = transSid(idx1);
            descrip = transDescrip(idx1);
            def = transDef(idx1);
            delG = transDelG(idx1);
            delGerror = transDelGerror(idx1);
            ec = transEC(idx1);
            if ~any(contains(presRxns, sid))% && rank ~= 1

                ub = transUB(idx1);
                lb = transLB(idx1);
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
                                rct = [rct,(reactstr+comptag)];
                                Rrecon(idx2) = reactstr+comptag;
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
                                rct = [rct,(reactstr+comptag)];
                                Lrecon(idx2) = reactstr+comptag;
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
                        outRxns{cnt,6} = "(T)";
                        outRxns{cnt,7} = traceNum;
                        outRxns{cnt,8} = "None";
                        outRxns{cnt,11} = descrip;
                        outRxns{cnt,12} = def;
                        outRxns{cnt,13} = delG;
                        outRxns{cnt,14} = delGerror;
                        outRxns{cnt,15} = ec;
                        outRxns{cnt,16} = "N/A";

                end
            end
        end
    end
end
end


