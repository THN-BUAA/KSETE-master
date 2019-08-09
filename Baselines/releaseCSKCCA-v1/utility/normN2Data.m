function [train_data,train_label,test_data,test_label] = normN2Data(source, target, idx)
 
ds = size(source,1); % fec*num
dt = size(target,1);

% sort data by each class
slabel = source(ds,:);
slabel(slabel>1) = 1;
posind = find(slabel == 1);
negind = find(slabel == 0);
temp1 = source(:,posind);
temp2 = source(:,negind);
source = [temp1, temp2];

% get training data, validation data
ratio = 0.9; % select 90% data of training data
trIdxPos = idx(1:floor(ratio*length(posind)));
valIdxPos = setdiff(idx(1:length(posind)),trIdxPos);
trIdxNeg = idx(length(posind)+1:length(posind)+floor(ratio*length(negind)));
valIdxNeg = setdiff(idx(length(posind)+1:end),trIdxNeg);

trIdx = [trIdxPos,trIdxNeg];
train_data = source(1:ds-1,trIdx);
train_label = source(ds,trIdx);
train_label(train_label>1) = 1;

% valIdx = [valIdxPos,valIdxNeg];
% val_data = source(1:ds-1,valIdx);
% val_label = source(ds,valIdx);
% val_label(val_label>1) = 1;
clear posind negind temp1 temp2

% get test data
tlabel = target(dt,:);
tlabel(tlabel>1) = 1;
posind = tlabel == 1;
negind = tlabel == 0;
temp1 = target(:,posind);
temp2 = target(:,negind);
target = [temp1, temp2];
clear posind negind temp1 temp2

test_data = target(1:dt-1,:);
test_label = target(dt,:);
test_label(test_label>1) = 1;

% N2 normalization
train_data = zscore(train_data,0,2);
% val_data = zscore(val_data,0,2);
test_data = zscore(test_data,0,2);   

end
