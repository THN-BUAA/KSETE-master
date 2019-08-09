clear all
close all
clc

addpath('.\utility\');
addpath('.\liblinear\'); % LR classifier

% you can change the path
load('.data\randmark_source_30_same.mat');
path = '.\data\HDPProject\';
savePath = '.\result_CSKCCA_dim70_Rep20\';

filename = dir(path);
filelen = length(filename);

for i = 1:filelen-2  % dataset
    subPath = [path,filename(i+2).name];
    if isdir(subPath)
        fn = dir(subPath);
        fnlen = length(fn);

        Rep = 20; % repeat
        for loop = 1:Rep
            sp = [savePath,filename(i+2).name];
            if exist(sp,'dir') == 0
                mkdir(sp);
            end

            score = []; prelabel = []; str = [];
            for j = 1:fnlen-2  % prdiction case 
                str{j,1} = fn(j+2).name;
                load([subPath, '\', str{j,1}]);
                sn = [sp, '\', str{j,1}];

                rdm = randmark{i,1}{j,:};
                [train_data,train_label,test_data,test_label] = normN2Data(source,target,rdm(loop,:));
                [score(j,:),prelabel(j,:)] = CSKCCA(train_data,train_label,test_data,test_label); 
            end

            Result_ONE = getMeasures(h,test_label);
            save([sp,'\',num2str(loop), '.mat'], 'score','prelabel','Result_ONE','str');
        end
    end
end

disp('done !')


