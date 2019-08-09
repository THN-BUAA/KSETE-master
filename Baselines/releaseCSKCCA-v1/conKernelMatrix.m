function [G1,G2] = conKernelMatrix(X1,X2,kernel1,kernel2,lrank)

N = size(X1,1);	% number of data1
NN = size(X2,1); % number of data2 

% get incompletely decomposed kernel matrices. K1 \approx G1*G1'
G1 = km_kernel_icd(X1,kernel1{1},kernel1{2},lrank);
G2 = km_kernel_icd(X2,kernel2{1},kernel2{2},lrank);

if size(G1,2)>size(G2,2)
    G1 = G1(:,1:size(G2,2));
elseif size(G1,2)<size(G2,2)
    G2 = G2(:,1:size(G1,2));
end  

% remove mean. avoid standard calculation N0 = eye(N)-1/N*ones(N);
G1 = G1-repmat(mean(G1),N,1);
G2 = G2-repmat(mean(G2),NN,1); 

