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
        if ~any(contains(presRxns, sid))% && rank ~= 1
            kid = supKid(idx1);
            ub = supUB(idx1);
            lb = supLB(idx1);
            gene = supGene(idx1);
            if ub > 0 && lb < 0
                rev = 1;
            else
                rev = 0;
            end
            if contains(rxn,"<=>")
                spl = strip(strsplit(rxn,"<=>"));
                Lstr = spl(1);
                Rstr = spl(2);
                dir = " <=> ";
            elseif contains(rxn,"=>") || contains(rxn,"-> ")
                spl = strip(strsplit(rxn,"=>"));
                Lstr = spl(1);
                Rstr = spl(2);
                dir = " => ";
            elseif contains(rxn," <=") || contains(rxn," <-")
                %Not the case for utilized database, consider ammending if time
            else
                warning("Master Database has unacceptable formagt for reaction '" + sid + "'")
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
    %             outRxns{cnt,8} = listGPR;
    %             outRxns{cnt,9} = listOrtho;
    %             outRxns{cnt,10} = listSpecies;

            end
        end
    end
end







%%%Function End%%%
end


%all(arrayfun(@(x) any(contains(viableMets,x)),rct))  
%Formatting for gene, target species ortholog, and ortholog origin is as
%follows:
%--Each gene that has a corresponding target species' gene is given
%   contained within a set of square brackets (if it is a set of genes
%   that are required, the whole set is given in square brackets
%   example: [ortholog1][ortholog2][ortholog3]
%--Each target species gene that corresponds to the ortholog is given in a
%   corresponding set of square brackets and separated by an 'or'. If
%   multiple genes are required, each set of target species genes that
%   correspond to that ortholog are given in parentheses and separated
%   internally by 'or'
%   example: [gene1.1 or gene 1.2][(gene 2.1) and (gene 2.2 or gene2.3)][(gene3.1) and (gene3.2)]
%--Origin species for each ortholog gene is stored in a set of square
%   brackets
%   example: [species1][species2][species3]
%--Each set of square brackets in one list corresponds to the same set of
%   brackets in the other
%   example:
%   [ortholog1]                         [ortholog2]                 [ortholog3]
%       /\                                  /\                           /\
%       ||                                  ||                           ||
%       \/                                  \/                           \/
%   [gene1.1 or gene 1.2][(gene 2.1) and (gene 2.2 or gene2.3)][(gene3.1) and (gene3.2)]
%       /\                                  /\                           /\
%       ||                                  ||                           ||
%       \/                                  \/                           \/
%   [species1]                          [species2]                    [species3]