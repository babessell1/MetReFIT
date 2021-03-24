function [outRxns] = tracerMasterReactExtract(met,masterDB,traceNum, presRxns)
keggFlag = 0;
mastRxns = string(masterDB.equation);
mastSid = string(masterDB.id);
masterDir = string(masterDB.direction);
%masterKid = string(masterDB.Kegg);
outRxns = {};

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
            elseif contains(rxn,"<=") || contains(rxn,"<-")
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
            end
        end
    end
end



%%%Function End%%%
end

%all(arrayfun(@(x) any(contains(viableMets,x)),rct))  
