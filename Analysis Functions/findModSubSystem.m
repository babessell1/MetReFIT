function [mode1, mode2, mode3, mode4, mode5] = findModSubSystem(metSubSystems)

subSystems = cell(0,2);

for i = 1:length(metSubSystems)
    metSubSysList = metSubSystems(i).subSystem;
    for j = 1:length(metSubSysList)
        metSubSys = metSubSysList{j};
        if metSubSys ~= "" && metSubSys ~= "Not assigned" && ~any(strcmp(metSubSys,subSystems(:,1)))
            subSystems{end+1,1} = metSubSys;
            subSystems{end,2} = 1;
        elseif metSubSys ~= "" && metSubSys ~= "Not assigned"
            numIdx = find(strcmp(metSubSys,subSystems(:,1)));
            for k = 1:length(numIdx)
                cntIdx = numIdx(k);
               subSystems{cntIdx,2} = subSystems{cntIdx,2} + 1;
            end
        end 
    end
end

subSystemCnts = cell2mat(subSystems(:,2));

[~,cntIdx] = max(subSystemCnts(:));
mode1 = subSystems{cntIdx,1};
newIdx = (ones(1, length(subSystems(:,2))));
newIdx(cntIdx) = 0;
subSystems = subSystems(logical(newIdx),:);

subSystemCnts = cell2mat(subSystems(:,2));

[~,cntIdx] = max(subSystemCnts(:));
mode2 = subSystems{cntIdx,1};
newIdx = (ones(1, length(subSystems(:,2))));
newIdx(cntIdx) = 0;
subSystems = subSystems(logical(newIdx),:);

subSystemCnts = cell2mat(subSystems(:,2));

[~,cntIdx] = max(subSystemCnts(:));
mode3 = subSystems{cntIdx,1};
newIdx = (ones(1, length(subSystems(:,2))));
newIdx(cntIdx) = 0;
subSystems = subSystems(logical(newIdx),:);

subSystemCnts = cell2mat(subSystems(:,2));

[~,cntIdx] = max(subSystemCnts(:));
mode4 = subSystems{cntIdx,1};
newIdx = (ones(1, length(subSystems(:,2))));
newIdx(cntIdx) = 0;
subSystems = subSystems(logical(newIdx),:);

subSystemCnts = cell2mat(subSystems(:,2));

[~,cntIdx] = max(subSystemCnts(:));
mode5 = subSystems{cntIdx,1};


