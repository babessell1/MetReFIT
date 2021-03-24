function [addRxns] = addTransport(met, traceNum, viableMets, excludeRxns)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%Determine if transfer reaction is needed, transport from cytosol if
%possible, if found in one other organelle other than cytosol,
%transport to cytosol first, if found in more than one organelle and
%NOT cytosol,pass and find if identical reactions producing metabolite
%are occuring in each organelle, if yes, add reaction to organelle and
%cytosol
if contains(met,"_mr") 
    metspl = strsplit(met,"_");
    mettag = "_"+metspl(end);
    mcomp = extractBetween(mettag,"[","]"); %stores the compartment code ex. "_[c]"
    nakedMet = metspl(1:end-1); %naked met is metabolite id with the compartment code stripped off ex. "C00001"
    vmets = viableMets(find(contains(viableMets,nakedMet))); %finds viable metabolites in all compartments matching the one in question
    vstr = string(vmets);
    vspl = arrayfun(@(s) strsplit(s,"_"),vstr,'UniformOutput',false); %seperates metabolite ids from the compartment code
    vtag = [];
    vcomp = [];
    for i = 1:length(vspl)
        vt = vspl{i}(2);
        vtag = [vtag;vt];
        vcomp = [vcomp;extractBetween(vt,"[","]")]; %removes formating from compartment code ex. "_[c]" --> "c"
    end
else
    metspl = strsplit(met,{'[',']'});
    mettag = "["+metspl{2}+"]";
    mcomp = metspl{2}; %stores the compartment code ex. "_[c]"
    nakedMet = metspl{1}; %naked met is metabolite id with the compartment code stripped off ex. "C00001"
    vmets = viableMets(find(contains(viableMets,nakedMet))); %finds viable metabolites in all compartments matching the one in question
    vstr = string(vmets);
    vspl = arrayfun(@(s) strsplit(s,{'[',']'}),vstr,'UniformOutput',false); %seperates metabolite ids from the compartment code
    vtag = [];
    vcomp = [];
    for i = 1:length(vspl)
        vt = "[" + vspl{i}(2) + "]";
        vtag = [vtag;vt];
        vcomp = [vcomp;vspl{i}(2)]; %removes formating from compartment code ex. "_[c]" --> "c"
    end
end

if ~isempty(vcomp)
    vcnt = sum(~contains(vcomp,"c"));
    vcompExl = vcomp(~contains(vcomp,"c"));
    vtagExl = vtag(~contains(vtag,"[c]"));
    if any(contains(vcomp,"c")) && mcomp ~= "c" 
        cytotag = strrep(mettag,"["+mcomp+"]","[c]");
        transID = "Transport_c"+mcomp+"_"+met;
        if ~ismember(transID,excludeRxns)
            addRxns = {transID,nakedMet+cytotag,met,nakedMet+cytotag+" <=> "+met,1,"(T)",traceNum,"TRANSPORT", 1.0};
        end
    elseif mcomp == "c" && vcnt == 1
        transtag = strrep(mettag,"[c]","["+vcompExl+"]");
        transID = "Transport_c"+vcompExl+"_"+nakedMet+transtag;
        if ~ismember(transID,excludeRxns)
            addRxns = {transID,nakedMet+transtag,met,met+" <=> "+nakedMet+transtag,1,"(T)",traceNum,"TRANSPORT", 1.0};
        end
    elseif mcomp ~= "c" && vcnt == 1
        cytotag = strrep(mettag,"["+mcomp+"]","[c]");
        transtag = strrep(mettag,"["+mcomp+"]","["+vcompExl+"]");
        transID1 = "Transport_c"+vcompExl+"_"+nakedMet+transtag;
        outRxn1 = {transID1,nakedMet+transtag,nakedMet+cytotag,nakedMet+cytotag+" <=> "+nakedMet+transtag,1,"(T)",traceNum,"TRANSPORT", 1.0};
        transID2 = "Transport_c"+mcomp+"_"+met;
        outRxn2 = {transID2,nakedMet+cytotag,met,nakedMet+cytotag+" <=> "+met,1,"(T)",traceNum,"TRANSPORT", 1.0};
        if ~ismember(transID2,excludeRxns) && ~ismember(transID1,excludeRxns)
            addRxns = vertcat(outRxn1,outRxn2);
        end
    else
        addRxns = {};
    end
else
    addRxns = {};
end
end

