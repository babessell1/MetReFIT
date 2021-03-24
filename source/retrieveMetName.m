function [metName] = retrieveMetName(metID, masterMets, masterMetsOld, masterMetsID, masterMetsIDold)
    %Searching to see if metabolite ID is a new ID and contained in the
    %master database.
    disp(metID)

    %disp(masterMetsID)
    matchIdx = find(contains(masterMetsID,metID));
    %If no new ID is found, see if metabolite ID is an old ID contained in
    %the master database.
    if isempty(matchIdx)
        matchIdx = find(contains(masterMetsIDold,metID));
    end
    
    if ~isempty(matchIdx)
        if ~ismissing(masterMets(matchIdx))
            comp = extractBetween(metID, "[","]");
            metName = masterMets(matchIdx) + "[" + comp + "]";
        elseif ~isempty(masterMetsOld(matchIdx))
            comp = extractBetween(metID, "[","]");
            metName = masterMetsOld(matchIdx) + "[" + comp + "]";
        else
            metName = metID;
        end
    else
        metName = metID;
    end
end