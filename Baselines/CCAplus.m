function [ PD,PF,Precision,F1,AUC,Accuracy,G_measure,MCC,Popt20,IFA] = CCAplus( source, target, LOC, randSeed)
%CCAPLUS Summary of this function goes here: Implement CCA+ algorithm.
%   Detailed explanation goes here
% INPUTS:
%   (1) source - a n1*(d1+1) matrix, the last column is the label where 0 is the majority class and 1 denotes the minority class.
%   (2) target - a n2*(d2+1) matrix.
%   (3) LOC    - the number of lines in each module. 
%   (4) randSeed - the seed of random.
% OUTPUTS:
%   
%
% Reference: X. Jing, F. Wu, X. Dong, F. Qi, and B. Xu, "Heterogeneous
%       cross-company defect prediction by unified metric representation
%       and cca-based transfer learning" in Proceedings of the 2015 10th
%       Joint Meeting on Foundations of Software Engineering, ser. ESEC/FSE
%       2015. New York, NY, USA: ACM, 2015, pp. 496¨C507.

warning('off');

% Default settings
if ~exist('randSeed','var')||isempty(randSeed)
    randSeed = 0; % 
end

rand('seed', randSeed);

%% Z-score normalization
[temp,mu,std]=zscore(source(:,1:end-1));
source(:,1:end-1) = temp;
[temp,mu,std]=zscore(target(:,1:end-1));
target(:,1:end-1) = temp;
    
%% UMR
numSrcIns = size(source,1);
numSrcFea = size(source,2)-1;
numTarIns = size(target,1);
numTarFea = size(target,2)-1;
if size(source,2)==63 % AEEEM is source
    if size(target,2)==21 % PROMISE
        commonIdx = [18, 29, 16, 41, 35, 33]; % 6 common CK metrics: wmc,dit,noc,cbo,rfc,lcom
        source = [source(:,commonIdx), source(:,setdiff(1:numSrcFea, commonIdx)), zeros(numSrcIns, numTarFea-length(commonIdx)), source(:,end)];
        target = [target(:,1:length(commonIdx)), zeros(numTarIns, numSrcFea-length(commonIdx)), target(:,length(commonIdx)+1:end)];
    else % JIRA
        source = [source(:,1:(end-1)), zeros(numSrcIns, numTarFea),source(:,end)];
        target = [zeros(numTarIns, numSrcFea), target];  
    end
elseif size(source,2)==21 % PROMISE is source
    if size(target,2)==63 % AEEEM
        commonIdx = [18, 29, 16, 41, 35, 33]; % 6 common CK metrics: wmc,dit,noc,cbo,rfc,lcom
        source = [source(:,1:(end-1)), zeros(numSrcIns, numTarFea-length(commonIdx)), source(:,end)];
        target = [target(:,commonIdx), zeros(numTarIns, numSrcFea-length(commonIdx)), target(:,setdiff(1:(numTarFea+1), commonIdx))];
    else % JIRA
        source = [source(:,1:(end-1)), zeros(numSrcIns, numTarFea),source(:,end)];
        target = [zeros(numTarIns, numSrcFea), target]; 
    end
else % JIRA is source
    source = [source(:,1:(end-1)), zeros(numSrcIns, numTarFea),source(:,end)];
    target = [zeros(numTarIns, numSrcFea), target]; 
end

%% Resampling - CCA requires that two datasets have the same number of instances 
ns = size(source,1);
nt = size(target,1);
if ns < nt
    idx = randperm(nt, nt);
    target = target(idx, :); % shuffle
    idxSel = randperm(nt, ns);
    targetRS = target(idxSel,:);
    sourceRS = source;
elseif ns > nt
    idx = randperm(ns, ns);
    source = source(idx, :); % shuffle
    idxSel = randperm(ns, nt);
    sourceRS = source(idxSel,:);
    targetRS = target;
end

%% CCA
[Ps,Pt,r] = canoncorr(sourceRS(:,1:end-1), targetRS(:,1:end-1)); % MATLAB build-in function

%% Projected data
srcX = source(:,1:end-1);
srcY = source(:,end);
tarX = target(:,1:end-1);
tarY = target(:,end);
newSource = [srcX*Ps(:,1:10),srcY];
newTarget = [tarX*Pt(:,1:10),tarY];

%% KNN - training and prediction
knn = fitcknn(newSource(:,1:end-1), newSource(:,end));
[preLabel, prob, ~] = predict(knn, newTarget(:,1:end-1)); %
probPos = prob(:,2);


try
    [ PD,PF,Precision,F1,AUC,Accuracy,G_measure,MCC,Popt20,IFA] = Performance(newTarget(:,end), probPos, LOC); % Call self-defined Performance()
catch
    PD=nan;PF=nan;Precision=nan;F1=nan;AUC=nan;Accuracy=nan;G_measure=nan;MCC=nan;Popt20=nan;IFA=nan;
 end

