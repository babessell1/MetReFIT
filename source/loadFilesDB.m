%%%This function is designed to facilitate the process of database input
%%%for the user in a step-by-step process rather than through a set of
%%%instructions. User can also skip this to save and write in the tracer2
%%%function itself to save the files.
function [presentDB, superDB, orthoDB, masterDB, masterMetsDB] = loadFilesDB()
    %Providing introductory output to explain input formatting
    fprintf('Please input the names/pathways for each of these databases:\n\t - The Present Database(most likely the same reaction list used in the model \n\t')
    fprintf(' - A \"Super\" database of reactions in closely-related species\n\t')
    fprintf(' - An ortholog database containing orthologs between the model species and the closely-related species \n\t')
    fprintf(' - A master database containing all reactions from a broader range of related-species \n\t\t')
    fprintf('(i.e., If the model species is the plant maize, the master database might contain all known plant reactions)\n\n\t')
    fprintf('*For help on these databases see \"Database Creation\" in the help folder\n\n');
    
    escapeFlag = 0;
    errorFlag = 0;
    
    %Present database input and table creation
    fprintf('\nPresent Database file (NOTE: ensure that the two sheets are titled \"Metabolite List\" and \"Reaction List\":\n','s');
        [userPresentInput, errorFlag] = retrieveUserDB();
        if ~errorFlag
            while ~errorFlag && ~escapeFlag
                [presentDB, escapeFlag, errorFlag] =  xlsdb2tabledb(userPresentInput, 1);
            end
        end
    
    
    %Super database input and table creation    
    if ~errorFlag
        escapeFlag = 0;
        fprintf('\n\nSuper Database file:\n');
        [userSuperInput, errorFlag] = retrieveUserDB();
        if ~errorFlag
            while ~errorFlag && ~escapeFlag
                [superDB, escapeFlag, errorFlag] =  xlsdb2tabledb(userSuperInput, 0);
            end
        end
    end
    %Ortholog database input and table creation
    if ~errorFlag
        escapeFlag = 0;
        fprintf('\n\nOrtholog Database file:\n');
        [userOrthoInput, errorFlag] = retrieveUserDB();
        if ~errorFlag
            while ~errorFlag && ~escapeFlag
                [orthoDB, escapeFlag, errorFlag] =  xlsdb2tabledb(userOrthoInput, 0);
            end
        end
    end
    %Master database input and table creation
    if ~errorFlag
        escapeFlag = 0;
        fprintf('\n\nMaster Database file(NOTE: ensure that the two sheets are titled \"Metabolite List\" and \"Reaction List\":\n');
        [userMasterInput, errorFlag] = retrieveUserDB();
        if ~errorFlag
            while ~errorFlag && ~escapeFlag
                [masterDB, escapeFlag, errorFlag] =  xlsdb2tabledb(userMasterInput, 1);
                [masterMetsDB, escapeFlag, errorFlag] = xlsdb2tabledb(userMasterInput, 2);
            end
        end
    end
    
      %If error occurred, making sure to return all databases as empty so
    %function does not pull an error
    if errorFlag
        presentDB = {};
        superDB = {};
        orthoDB = {};
        masterDB = {};
        masterMetsDB = {};
    end
    
    function [userDBInput, errorFlag] = retrieveUserDB()
        userRedo = "";
        errorFlag = 0;
        while lower(userRedo) ~= "yes" && ~errorFlag
            userDBInput = input('\n--please enter pathway\n','s');
            %Confirming pathway entry
            fprintf('\nThis is the pathway you entered: \n%s\n\n', userDBInput)
            userRedo = input('Is this correct? (yes/no)\n','s');
            while lower(userRedo) ~= "yes" && lower(userRedo) ~= "no"
                userRedo = input('Is this correct? (yes/no)\n','s');
                if lower(userRedo) ~= "yes" && lower(userRedo) ~= "no"
                    fprintf("Error with your input--please enter yes or no\n\n");
                end
            end
            %Asking if user would like to input another file if file is not
            %correct
            if lower(userRedo) ==  "no" || exist(userDBInput,'file') ~= 2
                if lower(userRedo) ~= "no"
                    fprintf("\n***ERROR***The file you entered does not exist.\n")
                end
                userRedo = input('\nWould you like to input another file? (yes/no)\n','s');
                while lower(userRedo) ~= "yes" && lower(userRedo) ~= "no"
                    fprintf("\nError with your input--please enter yes or no\n\n");
                    userRedo = input('\nWould you like to input another file? (yes/no)\n','s');
                end
                if lower(userRedo) == "yes" %If yes, setting to ask for another file
                    userRedo = "";
                else                        %If not, setting to escape
                    errorFlag = 1;         
                end
            end
        end
    end

  

    function [userDB, escapeFlag, errorFlag] = xlsdb2tabledb(userDBInput, dbType)
        escapeFlag = 0;
        errorFlag = 0;
        %Creating table from database
            try 
                if dbType == 0    
                    userDB = readtable(userDBInput);
                elseif dbType == 1
                    userDB = readtable(userDBInput,'Sheet','Reaction List');
                end
                if dbType == 2
                    userDB = readtable(userDBInput,'Sheet','Metabolite List');
                end
                escapeFlag = 1;  %Set to escape
            catch %If fails to create table, asks if user would like to input another file, if not, sets to escape
                fprintf("There is an issue with the formatting of the file you entered, ensure you have selected the correct document and that it is formatted correctly and try again. \nIf the problem persists, please reference documentation on creating databases to ensure you are using the correct format.\n\n") 
                userErrorRepeat = "";
                while lower(userErrorRepeat) ~= "yes" && lower(userErrorRepeat) ~= "no" 
                    userErrorRepeat = input('\nWould you like to enter another file (yes/no)?\n','s');
                    if lower(userErrorRepeat) ~= "yes" && lower(userErrorRepeat) ~= "no"
                        fprintf("\nError with your input--please enter yes or no\n\n");
                    end
                end
                if lower(userErrorRepeat) == "no"
                    errorFlag = 1;
                end
            end
    end
end

