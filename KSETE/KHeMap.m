function [ Bs, Bt, BsFull, BtFull ] = KHeMap( SData, TDataX ,isLinear, isMultiKen)
%HEMAP Summary of this function goes here
% Detailed explanation goes here:
% INPUTS:
%   (1) SData - n1*(p+1) matrix, p denotes the number of features(metrics), the last column denotes class labels {0,1} where 0 represents non-defective.
%   (2) TDataX - n2*q matrix, q is number of features.
%   (3) isLinear - {0,1}, 1 - HeMap with linear transformation, 0 - HeMap with nonlinear transformation; by default isLinear = 0; 
% OUTPUTS:
%   (1) Bs - n*k matrix, k denotes the number of new features obtained by projecting.
%   (2) Bt - m*k matrix
%


%% Settings of hyper-parameters
k = 1; % Dimension of projected feature space
beta = 1;
theta = 0.5;

%%
if nargin == 2
    isLinear = 1; % By defult
end

%% Original number of instances in source and target datasets -- VERY important
oriNumS = size(SData,1);
oriNumT = size(TDataX,1);

%% Preprocess - Ensure two datasets have same number of instances by random sampling.
if size(SData,1)>size(TDataX,1)
    % idx = randperm(size(T,1),size(SData,1)-size(T,1));
	idx = randi(size(TDataX,1), size(SData,1)-size(TDataX,1), 1); % radni(MAX, M,N) returns an M-by-N matrix containing pseudorandom integer values drawn from the discrete uniform distribution on 1:IMAX.
    TDataX = [TDataX;TDataX(idx,:)];    
elseif size(SData,1)<size(TDataX,1)
    % idx = randperm(size(SData,1),size(T,1)-size(SData,1));
	idx = randi(size(SData,1), size(TDataX,1)-size(SData,1), 1); 
    SData = [SData;SData(idx,:)]; 
end

S = SData(:,1:end-1); % Remove class label

Ct = [];
Cs = [];

numKernel = 3;
sigma = [1.5, 2, 3];
gamma = [0.5, 1 1.5];
% gamma = 1;
%% Calculate Bt, Bs
if isLinear==0 % NONLINEAR HeMap
    
    %Step1.1: Calculate Ct
    Idx = kmeans(TDataX,2); % two clusters; Idx is a (size(T,1))*1 vector (the minimum is 1), cluster label of each instance.    
	
    Ct = zeros(size(TDataX,1),2);
    if sum(Idx==1)>sum(Idx==2)
        Ct(find(Idx==1),1) = 1; % majority class
        Ct(find(Idx==2),2) = 1;
    else
        Ct(find(Idx==1),2) = 1; % majority class
        Ct(find(Idx==2),1) = 1;
    end
    
	
% 	Ct = T(find(Idx==1),:);
% 	Ct = [Ct;T(find(Idx==2),:)];
	
    %Step1.2: Calculate Cs
    label = SData(:,end);
    idx0 = find(label==0); % non-defective 
	idx1 = find(label==1);
    
    Cs = zeros(size(S,1),2);
    Cs(idx0,1) = 1;
    Cs(idx1,2) = 1;
	
%     Cs = [S(idx0,:);S(idx1,:)]; 
    

%     %Step2: case1 - Construct A without kernel
%     A1 = 2*theta^2*T*T' + beta^2/2*S*S' + (1-theta)*(beta+2*theta)*Ct*Ct';
%     A2 = beta*theta*(T*T'+S*S');
%     A3 = A2;
%     A4 = 2*theta^2*S*S' + beta^2/2*T*T' + (1-theta)*(beta+2*theta)*Cs*Cs';
    
    %Step2: case2 - Construct A with kernel. 
    rbf = sum(var(S));
    rbf2 = sum(var(TDataX));
    %rbf = 10;
    Kss = exp(-dist(S') / rbf);
    Ktt = exp(-dist(TDataX') / rbf2);  
    A1 = 2*theta^2*Ktt + beta^2/2*Kss + (1-theta)*(beta+2*theta)*Ct*Ct';
    A2 = beta*theta*(Ktt+Kss);
    A3 = A2;
    A4 = 2*theta^2*Kss + beta^2/2*Ktt + (1-theta)*(beta+2*theta)*Cs*Cs';
    
    
    A = [A1,A2;A3,A4];
    clear A1 A2 A3 A4;
    
    % Obtain tri-diagonal matrix T by lanczos approach
    [Td,V] = lanczos(A);
    % [Td,V] = lanczos(A,k);
    
    % Perform SVD on T
    [W,S1,V1] = svd(Td); % T = W*S*V1', when T is a symmetric matrix, then W equals to V1. Eigenvalues (elements of main diagonal in S1) are given in descending order.
    
    % Eigenvectors of A
    U = V*W;
    
    % Calculate Bt,Bs
    Bt = U(1:size(S,1),1:k);
    Bs = U(size(S,1)+1:end,1:k);
    
else % 
    % Case1:kernal based
    rbf = sum(var(S));
    rbf2 = sum(var(TDataX));
    
    if isMultiKen % Multiple kernels   
        Kss = exp(-dist(S').^2 / (gamma(1)*rbf));
        Ktt = exp(-dist(TDataX').^2 / (gamma(1)*rbf2));
        for i=2:length(gamma)
            Kss = Kss + exp(-dist(S').^2 / (gamma(i)*rbf));
            Ktt = Ktt + exp(-dist(TDataX').^2 / (gamma(i)*rbf2));
        end
        Kss = Kss / length(gamma);
        Ktt = Ktt / length(gamma);
                
    else
        % Single kernel
        Kss = exp(-dist(S') / rbf);
        Ktt = exp(-dist(TDataX') / rbf2);
    end  
    
    
    
    A1 = 2*Ktt + (beta^2)/2*Kss;
    A2 = beta*(Kss+Ktt);
    A3 = A2;
    A4 = 2*Kss + (beta^2)/2*Ktt;
    
    clear Kss Ktt;
    
%     % Case2: without kernel
%     A1 = 2*T*T' + (beta^2)/2*S*S';
%     A2 = beta*(S*S'+T*T');
%     A3 = A2;
%     A4 = 2*S*S' + (beta^2)/2*T*T';



    A = [A1,A2;A3,A4]; % A is a symmetric matrix
    clear A1 A2 A3 A4;
    
    % Obtain tri-diagonal matrix T by lanczos approach
    [Td,V] = lanczos(A);
    % [Td,V] = lanczos(A,k);
    clear A;
    
    % Perform SVD on T
    [W,~,~] = svd(Td);
    % [W,S1,V1] = svd(Td); % T = W*S*V1', when T is a symmetric matrix, then W equals to V1. Eigenvalues (elements of main diagonal in S1) are given in descending order.
    
    % Eigenvectors of A
    U = V*W;
    
    % Calculate Bt,Bs
    Bt = U(1:size(S,1),1:k);
    Bs = U(size(S,1)+1:end,1:k);
    clear U V W;
    
end

%% Recover original number of samples - VERY important!!!
BsFull = Bs; BtFull = Bt;
Bs = Bs(1:oriNumS,:);
Bt = Bt(1:oriNumT,:);

end

