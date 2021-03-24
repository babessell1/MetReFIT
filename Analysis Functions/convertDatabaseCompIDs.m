function [newMasterDB_rxnList, newMasterDB_metList] = convertDatabaseCompIDs(filename, dbType)



    if dbType == "super"
%         Not implemented, super database has already eliminated '_mr' tag
    elseif dbType == "master"
        dbTable = readtable(filename, 'Sheet', 'Reaction List');
        for i = 1:length(dbTable.Abbreviation)
            newRxnID = "";
            newRxnOldID = "";
            newRxn = "";
            
            rxnID = string(dbTable.Abbreviation(i));
            rxnOldID = string(dbTable.OldAbbreviation(i));
            rxn = string(dbTable.Reaction(i));
            
            splitRxnID = strsplit(rxnID, '_mr');
            for j = 1:length(splitRxnID)
                newRxnID = newRxnID + splitRxnID{j};
            end
            splitRxnOldID = strsplit(rxnOldID, '_mr');
            for j = 1:length(splitRxnOldID)
                newRxnOldID = newRxnOldID + splitRxnOldID{j};
            end
            splitRxn = strsplit(rxn, '_mr');
            for j = 1:length(splitRxn)
                newRxn = newRxn + splitRxn{j};
            end
            dbTable.Abbreviation{i} = newRxnID;
            dbTable.OldAbbreviation{i} = newRxnOldID;
            dbTable.Reaction{i} = newRxn;
        end
        newMasterDB_rxnList = dbTable;
        dbTable = readtable(filename, 'Sheet', 'Metabolite List');
        for i = 1:length(dbTable.Abbreviation)
            newMetID = "";
            newMetOldID = "";
            
            
            metID = string(dbTable.Abbreviation(i));
            metOldID = string(dbTable.old_abbrev(i));
            
            splitMetID = strsplit(metID, '_mr');
            for j = 1:length(splitMetID)
                newMetID = newMetID + splitMetID{j};
            end
            splitMetOldID = strsplit(metOldID, '_mr');
            for j = 1:length(splitMetOldID)
                newMetOldID = newMetOldID + splitMetOldID{j};
            end
            
            dbTable.Abbreviation{i} = newMetID;
            dbTable.old_abbrev{i} = newMetOldID;

        end
        newMasterDB_metList = dbTable;
    elseif dbType == "present"
    end


end