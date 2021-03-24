function [addRxns] = addFromUniversal(met, traceNum, viableMets, excludeRxns, masterDB, presentDB)

addRxns = {}; %list of reactions whose products contain the metabolite (met) in question
presRxns = string(presentDB.Abbreviation);

%Check for reactions in supermodel of which all reactants are viable
[uRxns] = tracerMasterReactExtract(met, masterDB, traceNum, presRxns); %pulls out reactions in the constructed "super" C4 plant model that produce the met in question
sprt = []; %products
if ~isempty(uRxns)
    for i = 1:length(uRxns(:,1))
        try
            rt = uRxns{i,2}; %reactants
            pt = uRxns{i,3}; %products
            rev = uRxns{i,5}; %reversibility
        catch
        end
        rx = uRxns{i,1}; %reaction ids
        sprt = [sprt,uRxns{i,2}];
        if all(ismember(rt,viableMets)) && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) %passes if all reactants of the reaction are viable and the reaction is not on the list of reactions to exclude
            uRxns{i, 9} = 0;
            blockPdt = sum(ismember(pt,blockedMets)); 
            uRxns{i, 10} = viaPdt;
            uRxns{i,17} = "None";
            addRxns = vertcat(addRxns,uRxns(i,:));
            
        elseif all(ismember(pt,viableMets)) && ~any(ismember(rx,excludeRxns)) && ~any(ismember(rx,presRxns)) && rev == 1
            uRxns{i, 9} = 0;  
            blockPdt = sum(ismember(pt,blockedMets));
            uRxns{i, 10} = viaPdt;
            uRxns{i,17} = "None";
            addRxns = vertcat(addRxns,uRxns(i,:));       
        end
       
    end
end

try
    addRxns = sortrows(addRxns, 10, 'descend');
catch
end
end

function [outRxns] = tracerMasterReactExtract(met,masterDB,traceNum, presRxns)
keggFlag = 0;
mastRxns = string(masterDB.equation);
mastSid = string(masterDB.id);
mastDir = string(masterDB.direction);
mastDescrip = string(masterDB.name);
mastEC = string(masterDB.ec_number);
mastDelG = masterDB.deltag;
mastDelGerror = masterDB.deltagerr;
mastPath = string(masterDB.pathways);
%masterKid = string(masterDB.Kegg);
outRxns = {};

metspl = strsplit(met,"[");
nakedMet = metspl(1);
tag = "[" + extractBetween(met,"[","]") + "]";
%determine rxns producing met of interest
idxE = find(contains(mastRxns,nakedMet));
if ~isempty(idxE)
    cnt = 0;
    for i = 1:length(idxE)
        idx1 = idxE(i);
        clear rxnRecon;
        pdt = [];
        rct = [];
        rxn = mastRxns(idx1);
        sid = mastSid(idx1);   
        descrip = mastDescrip(idx1);
        delG = mastDelG(idx1);
        delGerror = mastDelGerror(idx1);
        ec = mastEC(idx1);
        path = mastPath(idx1);
        
        
        if ~any(contains(presRxns, sid))
            if masterDir(idx1) == "="
                ub = 1000;
                lb = -1000;
            else
                ub = 1000;
                lb = 0;
            end
            if ub > 0 && lb < 0
                rev = 1;
            else
                rev = 0;
            end
            if contains(rxn," <=> ")
                spl = strip(strsplit(rxn,"<=>"));
                Lstr = spl(1);
                Rstr = spl(2);
                dir = " <=> ";
            elseif contains(rxn," => ") || contains(rxn," -> ")
                spl = strip(strsplit(rxn,"=>"));
                Lstr = spl(1);
                Rstr = spl(2);
                dir = " => ";
            elseif contains(rxn," <= ") || contains(rxn," <- ")
                spl = strip(strsplit(rxn,"<="));
                Rstr = spl(1);
                Lstr = spl(2);
                dir = " => ";
            else
                warning("Master Database has unacceptable format for reaction '" + sid + "'")
                continue
            end
            flipFlag=0;
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
            if exist('rxnRecon','var')
                cnt = cnt+1;                
                outRxns{cnt,1} = sid + tag;
                outRxns{cnt,2} = rct;
                outRxns{cnt,3} = pdt;
                outRxns{cnt,4} = strip(rxnRecon);
                outRxns{cnt,5} = rev;
                outRxns{cnt,6} = "(X)";
                outRxns{cnt,7} = traceNum;
                outRxns{cnt,8} = "N/A";
                outRxns{cnt,11} = descrip;
                outRxns{cnt,12} = "N/A";
                outRxns{cnt,13} = delG;
                outRxns{cnt,14} = delGerror;
                outRxns{cnt,15} = ec;
                outRxns{cnt,16} = path;              
            end
        end
    end
end

%%%Function End%%%
end