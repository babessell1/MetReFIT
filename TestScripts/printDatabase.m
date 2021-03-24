function printDatabase(dbVar, newFile)

userFile = fopen(newFile, 'w');

%Retrieving headers from table variable
varNames = dbVar.Properties.VariableNames;

for i = 1:length(varNames)
    if i == 1
        fprintf(userFile, '%s', varNames{i});
    else
        fprintf(userFile, '\t%s', varNames{i});
    end 
end
for i = 1:height(dbVar)
    for j = 1:length(varNames)
        if ~ismissing(dbVar(i,j))
            tableValue = string(dbVar{i,j});
            if j == 1
                fprintf(userFile, '\n%s', tableValue);
            else
                fprintf(userFile, '\t%s', tableValue);
            end
        else
            if j == 1
                fprintf(userFile, '\n');
            else
                fprintf(userFile, '\t');
            end
        end
    end
end

fclose(userFile);

