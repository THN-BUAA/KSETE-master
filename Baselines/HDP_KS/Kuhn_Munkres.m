function [M, MaxZjpp] = Kuhn_Munkres(A)
%KUHN_MUNKRES Summary of this function goes here
% Detailed explanation goes here: Use KM agorithm to complete the maximum weighted bipartile matching.
% INPUTS:
%   (1) A       - a n*n matrix, A(i,j) denotes a specific measure (i.e., weight) between i-th feature in
%   (X1,X2,...,Xn) and j-th feature in (Y1,Y2,...,Yn).
% OUTPUTS:
%   (1) M       - a n*n matrix, M(i,j)=1 means that Xi and Yj are matched,
%   M(i,j)=0 denotes Xi and Yj are not matched.
%   (2) MaxZjpp - the sum of weights for all matched pairs.
%



if size(A,1)~=size(A,2)
    error('For A, the number of rows must be equal to the number of columns!');
end

% clear all;
% global n;
% global nx;
% global ny;
global lx;
global ly;
global g;
global  slack;
global visx;
global visy;
% global INF;
global linker;

% g=[ 3 5 5 4 1;
%     2 2 0 2 2 ;
%     2 4 4 1 0	;
%     0 1 1 0 0;
%     1 2 1 3 3];

g = A;
INF = 7000000000;
n = size(A, 1);
nx = size(A, 1); 
ny = size(A, 2);
MAXLen = n; % 400

%%% KM algorithm O(n3)

for i = 1:MAXLen
    linker(i) = -1;
end

for i = 1:10000
    lx(i) = 0;
end
%lx[]=0;

for i = 1:10000
    ly(i) = 0;
end
for i = 1:nx
    
    lx(i) = -1 * inf;
    %  for(int j = 1; j <= ny; j++)
    for j = 1:ny
        if g(i,j) > lx(i)
            lx(i) = g(i,j);
        end
    end
end


for i =1:nx
    %(int x = 1; x <= nx; x++)
    for j =1:ny
        slack(j) = INF;
    end
    deadloop = 0;
    while deadloop < n   % 1 
        deadloop = deadloop +1;
        for j = 1:MAXLen  % *5
            visx(j) = 0;
        end
        for j = 1:MAXLen
            visy(j) = 0;
        end
        %if(DFS(x)),break;
        % disp(x);
        m = dfs1(i); % 0 - ; 1 - 
        if m || i >= nx
            break;     
        end
        d = INF;
        for j =1:ny              
            if visy(j) == 0 && d > slack(j)
                d = slack(j);
            end
        end
        for j = 1:nx
            if visx(j)
                lx(j) =lx(j) - d;
            end
        end
        for j = 1:ny
            if visy(j)
                ly(j) = ly(j) +d;
            else slack(j) = slack(j) -d;
            end
        end
        
    end
end

minDim = min(size(A,1), size(A,2)); 
M = zeros(minDim, minDim);
for i=1:n  
    % M(linker(i), i) = 1;
    if linker(i) ~= -1
        M(linker(i), i) = 1;
    end
end

MaxZjpp = 0;
for i =1 :nx
    if linker(i) ~= -1
        MaxZjpp =MaxZjpp +g(linker(i),i);
    end
end

% M       % Display the best match matrix {0,1} where 1 denotes match, 0 not match. 
% MaxZjpp % 

end


%%%%%%%%%%%%%%%%%%%%%%%%
function res = dfs1(x)

global lx;
global ly;
global visx;
global visy;
global g;
global linker;
global  slack;

visx(x) = 1;
for i =1:size(g, 2) %
    
    if visy(i)
        continue;
    end
    tmp = lx(x) + ly(i) - g(x,i);
    if tmp == 0   %(x,y)
        visy(i) = true;
        %back = linker(y) == -1 ;
        if linker(i) == -1
            back =1;
        else
            back = 0;
        end
        % if linker(y)
        if linker(i) > 0
            res = dfs1(linker(i));
        end
        %  end
        if back || res
            linker(i) = x;
            res = true;
            return
        end    
    else if slack(i) > tmp   
            slack(i) = tmp;
        end
    end
end
res = false;
end
