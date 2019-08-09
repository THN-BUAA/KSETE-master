function [PD, PF, Precision, F1, AUC, Accuracy, G_measure, MCC,Popt20,IFA] = HDP_KS(source, target, LOC)
%HDP_KS Summary of this function goes here
% Detailed explanation goes here:
% INPUTS:
%   (1) source - a n1*(d1+1) matrix, d1 denotes the number of features, the last column is the labels.  
%   (2) target - a n2*(d2+1) matrix, d2 denotes the number of features, the last column is the labels.
% OUTPUTS:
%   (1) PD, PF,..., MCC
%
% Reference: Nam J , Fu W , Kim S , et al. Heterogeneous Defect Prediction[J]. IEEE Transactions on Software Engineering, 2018:1-1.
%
%


%% Hyper-parameters
K = 10;             % the number of neighbors for relief
selectRatio = 0.15; % Select the top 15% metrics according to Nam J
threshold = 0.05;

%% Data and label
nMet1    = size(source, 2) - 1;
x_source = source(:, 1:nMet1);
y_source = source(:, nMet1+1);

nMet2    = size(target, 2) - 1;
x_target = target(:, 1:nMet2);
y_target = target(:, nMet2+1);


%% Feature selecton - relieff
[ranked1, ] = relieff(x_source, y_source, K); % ranked1 is the indices of columns in X ordered by attribute importance, meaning ranked1(1) is the index of the most important predictor.
% x_source = x_source(:, 1:floor(ranked1*selectRatio));

[ranked2, ] = relieff(x_target, y_target, K); 
% x_target = x_target(:, 1:floor(ranked2*selectRatio));

numCommFea = min(floor(size(x_source,2)*selectRatio), floor(size(x_target,2)*selectRatio));

% Update x_source, x_target
x_source = x_source(:, ranked1(1:numCommFea));
x_target = x_target(:, ranked2(1:numCommFea));

%% Metric matches - Maximum weighted bipartile matching 

% Construct metric matching matrix A
A = zeros(numCommFea, numCommFea);
for i=1:numCommFea % x_source
    for j=1:numCommFea % x_target
        [h, p] = kstest2(x_source(:, i), x_target(:, j)); % KSAnalyzer
        % A(i, j) = p;
        if p <= eps
            A(i, j) = 0;
        else
            A(i, j) = p;
        end
    end
end

% Find the best metric matches
if length(unique(A)) > 1
    [M, MaxZjpp] = Kuhn_Munkres(A);
    [idxRow, idxCol] = find(M==1);
    idx1 = []; idx2 = [];
    for i = 1:length(idxRow) % KSAnalyzer
        [h, p] = kstest2(x_source(:, idxRow(i)), x_target(:, idxCol(i))); 
        if p >= 0.05
            idx1 = [idx1, idxRow(i)];
            idx2 = [idx2, idxCol(i)];
        end
    end
    x_source = x_source(:, idx1);
    x_target = x_target(:, idx2);
else
    x_source = x_source(:, 1);
    x_target = x_target(:, 2);
end


%% Logistic regression classifier
% Training
model =glmfit(x_source, y_source, 'binomial', 'link', 'logit');

% Prediction
probPos = glmval(model, x_target, 'logit');
predLabel = double(probPos >= 0.5);

try
    [ PD,PF,Precision,F1,AUC,Accuracy,G_measure,MCC,Popt20,IFA] = Performance( y_target, probPos, LOC);
catch
    PD=nan;PF=nan;Precision=nan; F1=nan;AUC=nan;Accuracy=nan;G_measure=nan;MCC=nan;Popt20=nan;IFA=nan;
end

end

