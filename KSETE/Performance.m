function perf = Performance( actual_label, probPos, loc)
%PERFORMANCE Summary of this function goes here: Calculate various performance
%   Detailed explanation goes here
% INPUTS:
%   (1) actual_label  - The actual label, a column vetor, each row is an instance's class label {0,1}, 0 denotes nondefective, 1 denotes defective.
%   (2) probPos       - The probability of being predicted as postive class, which has the same size as actual_label.
%   (3) loc           - Change size or the lines of codes in each module when having no change size.
% OUTPUTS:
%   perf: a strcut including PF,PF,..,MCC,Popt,IFA.

if numel(unique(actual_label)) < 1
    error('Please make sure that the true label ''actual_label'' must has at least two different kinds of values.');   
end

if length(actual_label)~=length(probPos)
    error('The dimensions of actual labels and predicted labels are not identitical.');
end

if nargin < 3
    error('Some parameters are lost!');
end

if sum(probPos>1)>0
    error('Please check 2nd parameter!');
end

predicted_label = double(probPos>=0.5); %

cf=confusionmat(actual_label,predicted_label); % Matlab build-in function

if numel(cf)==1
    if unique(actual_label)==1 && unqiue(predicted_label)==1 % only positive sampels and all they are classified correctly
        TP = cf; TN = 0; FP = 0; FN = 0; % 
    else % only negative sampels and all they are classified correctly
        TP = 0; TN = cf; FP = 0; FN = 0;
    end
else
    TP=cf(2,2);
    TN=cf(1,1);
    FP=cf(1,2);
    FN=cf(2,1);
end

%% various performance
Accuracy = (TP+TN)/(FP+FN+TP+TN);
PD=TP/(TP+FN);
PF=FP/(FP+TN);
Precision=TP/(TP+FP);
F1=2*Precision*PD/(Precision+PD);
[X,Y,T,AUC]=perfcurve(actual_label, probPos, '1');% 
G_measure = (2*PD*(1-PF))/(PD+1-PF);
MCC = (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));

%% Popt20%
Popt = CalculatePopt([actual_label, probPos, loc]);


%% IFA (>=1)

% Step1:
A = [reshape(actual_label,numel(actual_label),1),reshape(predicted_label, numel(predicted_label),1),reshape(loc, numel(loc),1)];

A_predPos = A(A(:,2)==1,:);
A_predNeg = A(A(:,2)==0,:);

% Step2:
A_predPos_sort = sortrows(A_predPos,3); % ascending order
A_predNeg_sort = sortrows(A_predNeg,3); %

A = [A_predPos_sort;A_predNeg_sort];

% Step3:
idxs = find(A_predPos_sort(:,1)==1&A_predPos_sort(:,2)==1);
if ~isempty(idxs) % two case: A_predPos is empty or A_predPos is not empty but all samples in it are misclassified.
    IFA = idxs(1); % 
else
    IFA = sum(A(:,1)==0); % The worst case
end

perf.PD=PD; perf.PF=PF; perf.Precision=Precision; perf.F1=F1; perf.AUC=AUC; perf.Accuracy=Accuracy; perf.G_measure=G_measure; perf.MCC=MCC; perf.Popt=Popt; perf.IFA=IFA;
end

