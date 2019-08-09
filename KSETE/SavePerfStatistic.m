function SavePerfStatistic( perfs, names, targetCell, filePath, boolOutputMedian, numDecP)
%SAVEPERFSTATISTIC Summary of this function goes here: Save each performance as an Excel file.
%   Detailed explanation goes here
% INPUTS:
%   (1) perfs      - a 1*N cell array where N is the number of combinations (e.g., source->target). 
%   (2) names      - a 1*3 cell where names{1} is a cell of performance names, names{2} is a cell array of unique target names, and names{3} is a cell array of models
%   (3) targetCell - A 1*N cell array. 
%   (4) filePath (By default, user's destop)              - Specify the saving path of generated excel files.
%   (5) boolOuntputMedian (By default, false, i.e., mean) - Whether output the results in the way of 'median'. 
%   (6) numDecp (default -3)                              - a negative integer which denotes the number of decimal place.
%


%% Deault values

if ~exist('filePath', 'var')||isempty(filePath)
    filePath = [getenv('UserProfile'),'\Desktop'];
end
if ~exist('numDecP', 'var')||isempty(numDecP)
    numDecP = -3; % Number of decimal place
end
if ~exist('outputMedian', 'var')||isempty(boolOutputMedian)
    boolOutputMedian = false; % output in form of mean
end

%% Fuse experiements having same target by mean or median 
mytable = tabulate(targetCell); %
group = unique(targetCell,'stable'); % MUST 'stable'!!!
perfModels = cell(1,numel(group));
for i=1:numel(group) % each target data  
    len = cell2mat(mytable(i,2));% how many experiments have same target (the target name is group{i})
    temp = 0;
    for j=2:i
        temp = temp + cell2mat(mytable(j-1,2));
    end
    beginIdx = temp+1;       % the begin index of i-th element in group
    endIdx = beginIdx+len-1; % the end index of i-th element in group
    
    temp0 = cell(1, numel(perfs{1}));
    
    % performance
    for j=beginIdx:endIdx % each experiment having same target
        for k=1:numel(perfs{1}) % each model
            if boolOutputMedian
                temp0{k} = [temp0{k}; median(perfs{j}{k})]; % the median value of all runs in terms of each performance measure of k-th model on j-th experiment
            else
                temp0{k} = [temp0{k}; mean(perfs{j}{k})]; % perfModels{j}{k}): the average performance in terms of all measures of k-th model on j-th dataset for all runs.
            end
        end     
    end 
    
    perfModels{i} = temp0; % temp0: a 1*m cell (m is the number of models) where each element is a ki*d matrix
    
end


%% Initialization
numPerfs = size(perfModels{1}{1}, 2);
meanData = cell(1, numPerfs); % meanData={[],[],...,[]} - number of elements is # of performance measures, each element is a m*k matrix.
stdData = cell(1, numPerfs);  % {[],[],...,[]}
medianData = cell(1, numPerfs);

%% NAN -> zero
for i=1:numel(perfModels) % each dataset, perfModels={{[],[],..,[]},{[],[],...,[]},...,{[],[],...,[]}}
    for j=1:numel(perfModels{1,i}) % each model   
        temp = perfModels{1,i}{1,j};
        temp(isnan(temp)) = 0;
        perfModels{1,i}{1,j} = temp;
    end
end

%% Calculate mean and std
for i=1:numel(perfModels) % each dataset, perfModels={{[],[],..,[]},{[],[],...,[]},...,{[],[],...,[]}}
    for j=1:numel(perfModels{1,i}) % each model       
        temp0 = [];
        temp1 = [];
        model = perfModels{1,i}{1,j};    
        temp0(1,:) = nanmean(perfModels{1,i}{1,j},1); % 列均值
        if size(perfModels{1,i}{1,j},1)==1
            temp1(1,:) = zeros(1,size(perfModels{1,i}{1,j},2)); % 一个实数的方差等于0
        else
            temp1(1,:) = nanstd(perfModels{1,i}{1,j}); 
        end
        
        for k=1:numel(temp0) % each performance
            meanData{k}(i,j) = roundn(temp0(1,k),numDecP); % decimal place, and 'roundn' has no effect on INTEGER
            stdData{k}(i,j) = roundn(temp1(1,k),numDecP);  % decimal place 
        end
    end
end

%% Median values
for i=1:numel(perfModels) % each dataset, perfModels={{[],[],..,[]},{[],[],...,[]},...,{[],[],...,[]}}
    for j=1:numel(perfModels{1,i}) % each model       
        temp0 = [];
        temp0(1,:) = median(perfModels{1,i}{1,j},1);
        for k=1:numel(temp0) % each performance
            medianData{k}(i,j) = roundn(temp0(1,k),numDecP); % decimal place 
        end
    end
end

%% Whether replace 'meanData' with 'medianData' or not?
if boolOutputMedian
    meanData = medianData;
end

%% Save
AZ = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
for i=1:numel(names{1}) % each performance
    %% Save as EXCEL
    str = [filePath,'\',names{1}{i},'_Statistic','.xlsx'];
    
    % Sheet1: mean
    % Row/Column Title of Sheet1
    rowTitle = ['Data',names{2},'Average']; %!!!
    colTitle = names{3}; % names of performance measures
    xlswrite(str,rowTitle','Sheet1',['A1:A',num2str(length(rowTitle))]); % row title; 
    xlswrite(str,colTitle,'Sheet1',['B1:',AZ{numel(colTitle)+1},'1']); % column title
    
%     xlswrite(str,['Data',names{2}]','Sheet1',['A1:A',num2str(numel(names{2})+1)]); % row title; 
%     xlswrite(str,names{3},'Sheet1',['B1:',AZ{numel(names{3})+1},'1']); % column title
    
    % Mean ± std
    contentMS = [];
    if size(perfModels{1}{1,1},1)>1 % the case that perform each model multiple times on each dataset
        for j=1:size(meanData{i},1) % Each dataset
            for k=1:size(meanData{i},2) % Each model
                % contentMS{j,k} = num2str(meanData{i}(j,k));
                if boolOutputMedian
                    contentMS{j,k} = num2str(meanData{i}(j,k));
                else
                    contentMS{j,k} = [num2str(meanData{i}(j,k)), '±', num2str(stdData{i}(j,k))];
                end
            end
        end
    else 
        for j=1:size(meanData{i},1) % Each dataset
            for k=1:size(meanData{i},2) % Each model
                contentMS{j,k} = num2str(meanData{i}(j,k));
            end
        end
    end
    
    % Average value across all datasets
    for k=1:size(meanData{i},2) % Each model
        
        % Average
        if boolOutputMedian % median
            contentMS{j+1,k} = num2str(roundn(nanmean(meanData{i}(:,k),1),numDecP)); % average median
        else % mean
            contentMS{j+1,k} = [num2str(roundn(nanmean(meanData{i}(:,k),1),numDecP)), '±', num2str(roundn(nanstd(meanData{i}(:,k)),numDecP))];
        end 
       
    end
    
    xlswrite(str,contentMS,'Sheet1',['B2:',AZ{numel(names{3})+1}, num2str(length(rowTitle))]); % content
    close;
end



disp('Have saved successfully!');
end

