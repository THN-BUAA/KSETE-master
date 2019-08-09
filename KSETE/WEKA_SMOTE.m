function balData = WEKA_SMOTE(data, final_ratio, n_neighbor, rand_state)
%WEKA_SMOTE Implement SMOTE by using WEKA
% INPUTS:
%   (1) data                     - n_samples * (n_features+1) where the last column is the label;
%   (2) n_neighbor (default 5)   - the number of nearest neighbors;
%   (3) final_ratio (1, default) - the ideal ratio of the number of minority samples to that of majority samples;
%   (4) rand_state (0, default)  - random seed;
% OUTPUTS:
%   balData                      - [data; synthetic_minority]

if final_ratio<=(sum(data(:,end)==1)/sum(data(:,end)==0))
    balData = data; % Do not need resampling
    return;
end

%% Default value
if ~exist('final_ratio', 'var')||isempty(final_ratio)
    final_ratio = 1; % 
end
if ~exist('rand_state','var')||isempty(rand_state)
    rand_state = 1;
end
if ~exist('n_neighbor','var')||isempty(n_neighbor)
    n_neighbor = 5;
end


dataX = data(:,1:(end-1));
dataY = data(:,end);
label = cell(size(dataX,1),1);
for i=1:size(dataY,1)
    if dataY(i)==0
        label{i} = 'No';
    else
        label{i} = 'Yes';
    end
end

for i=1:size(dataX,2)
    feaNames{i} = ['x', num2str(i)];
end

feaNames{end+1} = 'Defect';
neg_size = sum(dataY==0);
pos_size = sum(dataY~=0);

%pert = round((neg_size-pos_size)/pos_size)*100;
pert = round((neg_size*final_ratio-pos_size)/pos_size)*100;
if pert<0 % 
    balData = data;
    return;
end
resampled = matlab2weka('data', feaNames, [num2cell(dataX) label]);
resampled.setClassIndex(resampled.numAttributes()-1);
% javaaddpath('SMOTE.jar');
smote = javaObject('weka.filters.supervised.instance.SMOTE');

smote.setOptions(['-K ', num2str(n_neighbor), ' -S ', num2str(rand_state)]);
smote.setPercentage(pert);
smote.setInputFormat(resampled);
try
    S = weka.filters.Filter.useFilter(resampled, smote); % WEKA_SMOTE is preferred.
    [balData,~] = weka2matlab(S,[]);
catch
    [mat, ~] = weka2matlab(resampled, []);
    minoSam = mat(mat(:,end)~=0,:);
    synSample = SMOTE( minoSam, size(minoSam,1), pert, n_neighbor, rand_state);
    balData = [mat; synSample];
end

end

