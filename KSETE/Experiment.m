
%% An example:

%javaaddpath('weka.jar'); 

dataNames ={{'Lucene', 'ant-1.7', 'ant-1.7', 'xalan-2.4'},...
            {'camel-1.4','camel-1.4','EQ', 'EQ'}};
targetCell = {'camel-1.4','camel-1.4','EQ', 'EQ'};

% Set saving path of emperiment results
filePath = strcat([getenv('UserProfile')  '\Desktop'], '\RQ1'); % failePath = 'D:\Experiments';
if ~exist(filePath,'dir') % 
    mkdir(filePath);
end

% Experiment Settings 
runtimes = 30;    % the number of runnings of prediction model
percent_tt = 0.9; % the percentage of training data in source data

perfNames = {'PD','PF','Precision','F1','AUC','Accuracy','G-Measure','MCC','Popt','IFA'};
modelNames = {'KSETE'}; % modelNames = {'KSETE','HDP-KS','CCAplus','CTKCCA'};

expNames = [];
for i=1:numel(dataNames{1})
    expNames{i} = [dataNames{1}{i},'_',dataNames{2}{i}];
end
sources = [];
targets = [];

perfs = cell(1,numel(dataNames{1}));

import weka.filters.*; 
import weka.*;

% dataPath = 'E:\Documents\Experiments\KSETE\Datasets-inOne\'; % MUST set!!!
dataPath = '..\Datasets-inOne\'; % MUST set!!!

for d = 1:numel(dataNames{1}) % Each dataset
    
    disp(['Data: ',num2str(d),' / ',num2str(numel(dataNames{1}))]);
       
    % Load source dataset    
    file1 = java.io.File([dataPath,dataNames{1}{d},'.arff']);  % create a Java File object (arff file is just a text file)
    loader = weka.core.converters.ArffLoader;  % create an ArffLoader object
    loader.setFile(file1);  % using ArffLoader to load data in file .arff
    insts = loader.getDataSet; % get an Instances object
    insts.setClassIndex(insts.numAttributes()-1); %  set the index of class label
    [sources,featureNames,targetNDX,stringVals,relationName] = weka2matlab(insts,[]); %{XXX,YYY}-->{0,1}
    sources = [sources(:, 1:end-1), double(sources(:, end)>0)]; % If defects(i) > 0, then defects(i) = 1, otherwise defects(i) = 0. 
    idx_loc_sor = 0;
    if ismember('ck_oo_numberOfLinesOfCode',featureNames) % AEEEM
        idx_loc_sor = 26;
    elseif ismember('loc',featureNames)                   % PROMISE
        idx_loc_sor = 11;
    else
        idx_loc_sor = 28;
    end
    
    % Load target dataset
    file2 = java.io.File([dataPath,dataNames{2}{d},'.arff']);  % create a Java File object (arff file is just a text file)
    loader = weka.core.converters.ArffLoader;  % create an ArffLoader object
    loader.setFile(file2);  % using ArffLoader to load data in file .arff
    insts = loader.getDataSet; % get an Instances object
    insts.setClassIndex(insts.numAttributes()-1); %  set the index of class label
    [targets,featureNames,targetNDX,stringVals,relationName] = weka2matlab(insts,[]); %{false,true}-->{0,1}
    targets = [targets(:, 1:end-1), double(targets(:, end)>0)];
    idx_loc_tar = 0;
    if ismember('ck_oo_numberOfLinesOfCode',featureNames) % AEEEM
        idx_loc_tar = 26;
    elseif ismember('loc',featureNames)                  % PROMISE
        idx_loc_tar = 11;
    else
        idx_loc_tar = 28;
    end   
    
    % Remove instances have zero LOC
    sources(sources(:,idx_loc_sor)==0,:) = []; % sources(sources(idx_loc_sor)==0,:) = []
    targets(targets(:,idx_loc_tar)==0,:) = [];
    
    % Remove duplicated instances
    sources = unique(sources,'rows','stable');
    targets = unique(targets,'rows','stable');
    
    % Remove instances having missing values
    [idx_r idx_c] = find(isnan(sources));
    sources(unique(idx_r),:) = [];
    [idx_r idx_c] = find(isnan(targets));
    targets(unique(idx_r),:) = [];
    

    % Shuffle the instances  
    if runtimes >= 1
        rand('state',0);
        sources = sources(randperm(size(sources,1),size(sources,1)),:); % Disrupt the order of rows, which is beneficial to the following random resampling.
        targets = targets(randperm(size(targets,1),size(targets,1)),:);
    end
    
    LOC = targets(:,idx_loc_tar); 
    
    
    % Initialization
    PD_ksete=zeros(runtimes,1);PF_ksete=zeros(runtimes,1);Precision_ksete=zeros(runtimes,1);F1_ksete=zeros(runtimes,1);AUC_ksete=zeros(runtimes,1);
    Accuracy_ksete=zeros(runtimes,1);G_measure_ksete=zeros(runtimes,1);MCC_ksete=zeros(runtimes,1);Popt_ksete=zeros(runtimes,1);IFA_ksete=zeros(runtimes,1);  
    
    sourcesCopy = sources;  
    for i=1:runtimes  
        disp(['runtimes:',num2str(i),' / ',num2str(runtimes)]);
        % Obtain the training data
        rand('seed',i);
        idx = randperm(size(sourcesCopy,1),round(percent_tt*size(sourcesCopy,1)));  %
        trainData = sourcesCopy(idx,:);    
        while numel(unique(trainData(:,end)))==1 % Avoid the label has only one kind.
            idx = randperm(size(sourcesCopy,1),round(percent_tt*size(sourcesCopy,1)));
            trainData = sourcesCopy(idx,:);
        end
        sources = trainData;
        
        
        %% The proposed KSETE
        disp('KSETE...');       
        source = sources;
        target = targets;
        
        perf = KSETE(source, target,LOC); 
        PD_ksete(i,:)=perf.PD; PF_ksete(i,:)=perf.PF; Precision_ksete(i,:)=perf.Precision; F1_ksete(i,:)=perf.F1; Accuracy_ksete(i,:)=perf.Accuracy;
        AUC_ksete(i,:)=perf.AUC;G_measure_ksete(i,:)=perf.G_measure;MCC_ksete(i,:)=perf.MCC;Popt_ksete(i,:)=perf.Popt;IFA_ksete(i,:)=perf.IFA;
        
    end % End of runs
    
    perfs{d} = {[PD_ksete,PF_ksete,Precision_ksete,F1_ksete,AUC_ksete,Accuracy_ksete,G_measure_ksete,MCC_ksete,Popt_ksete,IFA_ksete]};
    
    save([filePath,'\perfs.mat'],'perfs'); 
end
 




