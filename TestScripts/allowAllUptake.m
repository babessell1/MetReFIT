function model = allowAllUptake(model)
rxns2change = [];
for i = 1:length(model.rxns)
    rxnID = model.rxns(i);
    if contains(rxnID, 'EX_') == 1
        rxns2change = [rxns2change;rxnID];
    end
end
model = changeRxnBounds(model,rxns2change,1000,'u');

