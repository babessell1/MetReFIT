function [metSubSystems] = retrieveMetSubSystems(blockedMetInfo)
    metInfoMets = blockedMetInfo(:,1);
    metInfoGeneSubSys = blockedMetInfo(:,4);
for i = 1:length(metInfoMets)
    metSubSystems(i).met = metInfoMets{i};
    metGeneSubSys = metInfoGeneSubSys{i};
    metSubSys = metGeneSubSys(:,2);
    subSysList = "";
    for j = 1:length(metSubSys)
        subSys = metSubSys{j};
        if iscell(subSys)
            for k = 1:length(subSys)
                if ~contains(subSysList,subSys{k})
                    subSysList = [subSysList;subSys{k}];
                end
            end
        else
            if contains(subSysList, subSys)
            else
                subSysList = [subSysList;subSys];
            end
        end
    end
    metSubSystems(i).subSystem = subSysList;
end
    
    
end