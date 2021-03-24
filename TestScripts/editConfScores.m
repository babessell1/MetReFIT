function [originData] = editConfScores(originData, userOp, originID, species, score)
    if userOp == "edit"
        if isnumeric(originID) && originID + 4 <= length(originData)
            if species == "origin"
                originData(originID + 4).origin = "Super Database: " + score;
            elseif species == "score"
                originData(originID +4).confScore = score;
            else
                warning('This is not an editable field');
            end
        elseif originID == "GPR" && isnumeric(score)
            originData(1).confScore = score;
        elseif originID == "master" && isnumeric(score)
            originData(2).confScore = score;
        elseif originID == "transport" && isnumeric(score)
            originData(3).confScore = score;
        elseif originID == "current" && isnumeric(score)
            originData(4).confScore = score;
        
        elseif (originID == "GPR" || originID == "master" || originID == "transport" || originID == "current") && ~isnumeric(score)
            warning('The name of the origin cannot be edited.')
        else
            warning('This is not an editable confidence score');
        end
    elseif userOp == "add"
        if isnumeric(score)
            originData(length(originData)+1).origin = "Super Database: "+ species;
            originData(length(originData)).confScore = score;
        else
            warning('This reaction origin cannot be added, try entering a number for the confidence score.');
        end
    %Removes unwanted origins and confidence scores
    elseif userOp == "remove"
        if isnumeric(originID) && originID+4 <= length(originData)
            idx = ones(1,length(originData));
            idx(originID+4) = 0;
            originData = originData(logical(idx));
        else
            warning('Confidence score cannot be edited or does not exist.')
        end
    end
    
end