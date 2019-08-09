function PoptValue = CalculatePopt(data,firstPercent)
%CALCULATEPOPT Summary of this function goes ahere
%   Detailed explanation goes here
% INPUTS:
%   (1) data         - n*3 matrix, where 3 columns are actual defects(or actual label {0,1} where 1 is defective),
%   predicted defects(or predicted probability of being positive class), and LOC, respectively. 
%   (2) firstPercent - 0.2, By default.
% OUTPUTS:
%   PoptValue        - return the Popt value.
%


% Set default value
if nargin == 1
    firstPercent = 0.2;
end
if firstPercent<0||firstPercent>1
    error('Parameter value crosses the border!');
end

if min(data(:,2))<0
    error('The 2nd column of data cannot include negative value!');
end

% Calculate actual defect density (NOTE:)
defDensity = (data(:,1))./(data(:,3)); %
data = [data,defDensity]; % Add a column, i.e., the actual defect sensity

% Calculate area of optimal model, actual model, and worst model.
area_m = CalculateArea(data, firstPercent, 'm');
area_optimal = CalculateArea(data, firstPercent, 'opt');
area_worst   = CalculateArea(data, firstPercent, 'worst');


% Calculate Popt 
try
    PoptValue = 1 - (area_optimal-area_m)/(area_optimal-area_worst);

catch
    PoptValue = nan;
end

end

function area = CalculateArea(data, firstPercent, type)
%CALCULATEAREA Summary of this function goes ahere
%   Detailed explanation goes here
% INPUTS:
%   (1) data - n*4 matrix, where 4 columns are RelDefects(or RelLabel), PredDefects(or ProbPos), LOC, defectDensity, respectively.
%   (2) firstPercent - [0,1], e.g., 0.2 denotes the first 20% LOC.
%   (3) type - {'opt','worst','m'}
% OUTPUTS:
%   area - 'approximate' area under a curve.


%%
if size(data,2)~=4
    error('Not enough columns in data!');
end

if strcmp(type, 'opt')
    % Dinvid into two parts - positive and negative, according to actual defects or label
    negData = data(data(:,1)==0,:);
    posData = data(~data(:,1)==0,:);
    
    % Sorting by actual defect density in descending order and then by LOC in ascending order if having duplicated defect density values 
    posData = sortrows(posData, [-4, 3]);
    negData = sortrows(negData, [-4, 3]);
    
    % Calculate the cumulative sum of LOC and actual defects
    tempData = [posData; negData]; % Owing to 'opt'
    cumLOC = cumsum(tempData(:,3));
    cumDefects = cumsum(tempData(:,1));
    
    % Calculate proportion
    Xs = cumLOC/cumLOC(end);
    Ys = cumDefects/cumDefects(end);
  
elseif strcmp(type,'worst')
    % Dinvid into two parts - positive and negative, according to actual defects or label
    negData = data(data(:,1)==0,:);
    posData = data(~(data(:,1)==0),:);
    
    % Sorting with actual defect density in ascending order and then LOC in descending order if having duplicated defect density values 
    posData = sortrows(posData, [4, -3]);
    negData = sortrows(negData, [4, -3]);
    
    % Calculate the cumulative sum of LOC and actual defects
    tempData = [negData;posData]; % Owing to 'worst', means that all samples are misclassified
    cumLOC = cumsum(tempData(:,3));
    cumDefects = cumsum(tempData(:,1));
    
    % Calculate proportion
    Xs = cumLOC/cumLOC(end);
    Ys = cumDefects/cumDefects(end);
    
else % 'm'

    % Add a column - predcited defect density
    if max(data(:,2))<=1
        predLabel = double(data(:,2)>=0.5); % if not smaller than 0.5, then regard as defective class
        predDefDensity = (predLabel)./data(:,3);
    else
        predDefDensity = (data(:,2))./data(:,3);
    end
    data = [data, predDefDensity];%
    % Divide into two parts - positive/negtive with the predicted defect
    if max(data(:,2))<=1 % probability of being positive
        prePosData = data(data(:,2)>=0.5,:);    % be defective (i.e., positive) if positive probability is not smaller than 0.5  
        preNegData = data(~(data(:,2)>=0.5),:); % non-defective
    else % defect counts
        preNegData = data(data(:,2)==0,:);
        prePosData = data(data(:,2)~=0,:);
    end
    % Sorting by predicted defect density in descending order, for the replicated sampels, sort in ascending order by LOC 
    preNegData = sortrows(preNegData, [-5, 3]);
    prePosData = sortrows(prePosData, [-5, 3]);
    % Combination
    data = [prePosData; preNegData];
        

    % Calculate the cumulative sum of LOC
    cumLOC = cumsum(data(:,3));
    
    
    % Calculate proportion
    Xs = cumLOC/cumLOC(end);
%     Ys = cumDefects/cumDefects(end);
    if max(data(:,2))<=1 % classification {0,1}
        cumDefects = cumsum(data(:,1).*double(data(:,2)>=0.5));
        Ys = cumDefects/sum(data(:,1));
    else % regression
        % temp = data(:,2);
        % temp(temp>2*data(:,1)) = 0; % 
        % actual = data(:,1); 
        % temp(temp>actual) = actual(temp>actual) - (temp(temp>actual) - actual(temp>actual)); % 
        % cumDefects = cumsum(temp);
        % Ys = cumDefects/sum(data(:,1)); % 
    end
end

% Identify the index of element which is nearest with 'firstPercent'
idx = find(Xs<=firstPercent);
if isempty(idx) % 
    subArea = 0.5*Xs(1)*Ys(1); 
else
    if Xs(idx(end))==firstPercent
        subArea = zeros(idx(end),1);                                                    % Initialization
        subArea(1) = 0.5*Xs(1)*Ys(1);                                                 % area of a triangle
        if idx(end)>1
            subArea(2:idx(end)) = 0.5*(Ys(1:(idx(end)-1))+Ys(2:idx(end))).*abs(Xs(1:(idx(end)-1))-Xs(2:idx(end))); % area of each trapezoid i.e, area = 0.5*(a+b)*h;
        end
    else
        subArea = zeros(idx(end)+1,1);                                                    % Initialization
        subArea(1) = 0.5*Xs(1)*Ys(1);                                                 % area of a triangle
        if idx(end)>1
            subArea(2:idx(end)) = 0.5*(Ys(1:(idx(end)-1))+Ys(2:idx(end))).*abs(Xs(1:(idx(end)-1))-Xs(2:idx(end))); % area of each trapezoid i.e, area = 0.5*(a+b)*h;
        end
        subArea(idx(end)+1) = 0.5 * (Ys(idx(end))+Ys(idx(end))+(((Ys(idx(end)+1)-Ys(idx(end)))*(firstPercent-Xs(idx(end))))/(Xs(idx(end)+1)-Xs(idx(end))))) * (firstPercent-Xs(idx(end))); % the area of a small trapezoid (the last trapezoid)
    end
end


% Calculate the sum of areas of all small regions
area = sum(subArea);

end