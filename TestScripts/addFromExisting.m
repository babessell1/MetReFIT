function [addRxns] = addFromExisting(met, traceNum, viableMets, excludeRxns, presentDB)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

addRxns = {}; %list of reactions whose products contain the metabolite (met) in question
[eRxns] = tracerExistReactExtract(met,presentDB,traceNum); %pulls out reactions already in the model that produce the met in question
%check for reactions already in model that are producing the metabolite of
%interest, skip mets whose production should already be feasible
if ~isempty(eRxns) && length(eRxns(:,1))>1
    for i = 1:length(eRxns(:,1))
        rt = eRxns{i,2}; %reactants
        pt = eRxns{i,3}; %products
        rx = eRxns{i,1}; %reaction id
        if ~ismember(rx,excludeRxns) && ~all(ismember(pt,viableMets)) %prevents reversible reactions causing the function to get stuck in an infinite loop
            if length(rt) - sum(ismember(rt,viableMets)) == 1
                eRxns{i,8} = "Present";
                eRxns{i,9} = 1.0;
                eRxns{i,10} = 1;
                addRxns = vertcat(addRxns, eRxns(i,:));
            end
        end
    end
end
try
    addRxns = sortrows(addRxns,10);
catch
end
end

function [outRxns] = tracerExistReactExtract(met,presentDB,traceNum)
% Function for determining if reacitons that already exist in the model are
% producing the met of interest.
% arguments: metabolite, database of existing reactions, tracing number
% from trace function (for tracking how far downstream from original
keggFlag = 0;
presRxns = string(presentDB.Reaction);
presAbb = string(presentDB.Abbreviation);
try
    presLB = presentDB.LowerBound;
    presUB = presentDB.UpperBound;
catch
    presLB = presentDB.Lowerbound;
    presUB = presentDB.Upperbound;
end
outRxns = {};
%determine rxns producing met of interest
idxE = find(contains(presRxns,met));
if ~isempty(idxE)
    cnt = 0;
    for i = 1:length(idxE)
        idx1 = idxE(i);
        clear rxnRecon;
        pdt = [];
        rct = [];
        rxn = presRxns(idx1);
        abb = presAbb(idx1);
        ub = presUB(idx1);
        lb = presLB(idx1);
        if ub > 0 && lb < 0
            rev = 1;
            dir = " <=> ";
        else
            rev = 0;
            dir = " => ";
        end
        if contains(rxn,"<=>")
            spl = strip(strsplit(rxn,"<=>"));
            Lstr = spl(1);
            Rstr = spl(2);
            %dir = " <=> ";
        elseif contains(rxn,"=>") || contains(rxn,"-> ")
            spl = strip(strsplit(rxn,{'=>','->'}));
            Lstr = spl(1);
            Rstr = spl(2);
            %dir = " => ";
        elseif contains(rxn,"<=") || contains(rxn,"<-")
            %Not the case for utilized database, consider ammending if time
        else
            warning("Master Database has unacceptable format for reaction '" + abb + "'")
            continue
        end
        flipFlag = 0;
        if lb < 0 && contains(Lstr,met)
            flipFlag = 1;
            prodList = strip(strsplit(Lstr," "));
            Lrecon = strings(size(prodList));
            for idx2 = 1:length(prodList)
                prodstr = prodList(idx2);
                if isnan(str2double(prodstr)) && ~startsWith(prodstr,"+") && ~startsWith(prodstr,"=") && ~startsWith(prodstr,"<") && ~startsWith(prodstr,"-")
                    if keggFlag
                        %%% space for reading Kegg format (not implemented)
                    else
                        pdt = [pdt,prodstr];
                        Lrecon(idx2) = prodstr;
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
                        rct = [rct,(reactstr)];
                        Rrecon(idx2) = reactstr;
                    end
                    
                else
                    Rrecon(idx2) = reactstr;
                end
            end
        elseif ub > 0 && contains(Rstr,met)
            prodList = strip(strsplit(Rstr," "));
            Rrecon = strings(size(prodList));
            for idx2 = 1:length(prodList)
                prodstr = prodList(idx2);
                if isnan(str2double(prodstr)) && ~startsWith(prodstr,"+") && ~startsWith(prodstr,"=") && ~startsWith(prodstr,"<") && ~startsWith(prodstr,"-")
                    if keggFlag
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                        pdt = [pdt,(prodstr)];
                        Rrecon(idx2) = prodstr;
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
                        rct = [rct,(reactstr)];
                        Lrecon(idx2) = reactstr;
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
        elseif isempty(Lrecon) && isempty(Rrecon)
        else
            error("Check tracerExistReactExtract.m")
        end
        if exist('rxnRecon','var')
            cnt = cnt+1;
            outRxns{cnt,1} = abb;
            outRxns{cnt,2} = rct;
            outRxns{cnt,3} = pdt;
            outRxns{cnt,4} = strip(rxnRecon);
            outRxns{cnt,5} = rev;
            outRxns{cnt,6} = "(0)";
            outRxns{cnt,7} = traceNum;
            outRxns{cnt,8} = "NONE";
            outRxns{cnt,9} = 1.0;
            
        end
    end
end
end

