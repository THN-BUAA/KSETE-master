function perf = KSETE( source, target,LOC, K, isMultiKen)  
% function [PD, PF, Precision, F1, AUC, Accuracy, G_measure, MCC, Popt, IFA] = KHeMapE( source, target,LOC, K, isMultiKen)  
%KHEMAPE Summary of this function goes here
%   Detailed explanation goes here: kernel HeMap ensemble 
% INPUTS:
%   (1) source - n1*(d1+1) matrix, each row is an instance, the last column is the label.
%   (2) target - n2*(d2+1) matrix
%   (3) LOC    - the lines of code of software modules of target dataset
%	(4) K      - Number of the base learners
%	(5) isMultiKen - 
% OUTPUTS:
%   perf - a struct of performance measures.
% 
% 


%% Hyper-parameters
if nargin==3
    K = 5; % Number of the base learners
    isMultiKen = true; % multiplr kernel by default
end

numSubFeat = size(source,2)-1;

%% Perform SMOTE
source = WEKA_SMOTE(source, 0.85, 5, 1); % 

%% Training
KHMcells = cell(1, K);
for i=1:K
    % disp(['i=', num2str(i)])
    %% subset of source
    idx1 = randi(size(source, 1), size(source, 1), 1); % % radni(MAX, M,N) returns an M-by-N matrix containing pseudorandom integer values drawn from the discrete uniform distribution on 1:IMAX.
    sub_source = source(idx1, :);
    
    %% subset of features
    idx2 = randi(size(source, 2)-1, 1, numSubFeat);
    sub_source = [sub_source(:, idx2), sub_source(:, end)];
    subModel.feature = idx2;
    %% Z-score
    [temp, mu, std] = zscore(sub_source(:,1:end-1));
    sub_source(:,1:end-1) = temp;
    [temp, mu, std] = zscore(target(:,1:end-1));
    target(:,1:end-1) = temp;
    
    %% Mapping
    [Bs, Bt, BsFull, BtFull] = KHeMap(sub_source, target(:,1:end-1), 1, isMultiKen); % 0 - nonlinear or kernel, 1 - linear or kernel
    
    source_new = [Bs, sub_source(:,end)];
    target_new = [Bt, target(:,end)];
    
    %% Logistic regression
    
    model =glmfit(source_new(:,1:end-1),source_new(:,end),'binomial', 'link', 'logit'); % 
    
    subModel.model = model;
    subModel.testData = target_new;
    subModel.FNorm = norm(BsFull-BtFull, 'fro');
    
    KHMcells{i} = subModel;
end

%% Testing
% case1: same weights
weights = ones(1, K) / K;


probPos = [];
for i=1:K
    model = KHMcells{i}.model;
    target_new = KHMcells{i}.testData;
    probPos(:, i) = glmval(model,target_new(:,1:end-1), 'logit'); % proPos denotes the predicted probability of being positive class   
end
finalProbPos = sum(bsxfun(@times, probPos, weights), 2);
predLabel = double(finalProbPos>=0.5);


try
    perf = Performance( target(:,end), finalProbPos, LOC);
catch
    perf.PD=nan; perf.PF=nan; perf.Precision=nan; perf.F1=nan; perf.AUC=nan; perf.Accuracy=nan; perf.G_measure=nan; perf.MCC=nan; perf.Popt=nan; perf.IFA=nan;
end
end

