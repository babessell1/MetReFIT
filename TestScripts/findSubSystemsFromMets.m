function [blockedMetInfo] = findSubSystemsFromMets(model, metList)
%Finds the associated subsystems for each blocked metabolite
%Work flow:
%   1) Finds reactions associated with each metabolite
%   2) Finds genes associated with these reactions
%   3) Finds subsystems for genes

blockedMetInfo = cell(length(metList), 4);

for i = 1:length(metList)
    geneList = [];
    subSystemList = [];
    blockedMetInfo{i,1} = metList(i);
    %Find reactions associated with the blocked metabolite
    [rxnList, rxnFormulaList] = findRxnsFromMets(model,metList(i));
    %Find genes associated with each rxn -- goes reaction by reaction so as to be used for an input for the next part
    for j = 1:length(rxnList)
        rxnGeneSysInfo(j).reaction = rxnList(j);
        rxnGeneSysInfo(j).rxnFormula = rxnFormulaList(j);
        
        [rxnGeneList] = findGenesFromRxns(model,rxnList(j));
        rxnGeneList = rxnGeneList{1};
        
        rxnGeneSysInfo(j).genes = rxnGeneList;
        
        geneList = vertcat(geneList, rxnGeneList);
        
        %Find subsystems associated with genes
        [geneSubSystems, singleList] = findSubSystemsFromGenes(model,rxnGeneList);
        
        rxnGeneSysInfo(j).subSystem = singleList;
        subSystemList = vertcat(subSystemList, singleList);
    end

    %Saving information in cell array
    blockedMetInfo{i,2} = rxnGeneSysInfo;
    blockedMetInfo{i,3} = geneList;
    blockedMetInfo{i,4} = subSystemList;
  
   
end