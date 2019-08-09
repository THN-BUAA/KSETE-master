function [score,predict_label,dim,s] = CSKCCA(train_data, train_label, test_data, test_label)
% cost-sentitive transfer KCCA

lrank = 70;  % number of components in incomplete Cholesky decompose
reg = 1E-5; % regularization

d1 = pdist(train_data');
sigma1 = mean(d1);
d2 = pdist(test_data');
sigma2 = mean(d2);

kernel1 = {'gauss',1/sigma1};   % kernel type and kernel parameter for data set 1
kernel2 = {'gauss',1/sigma2};   % kernel type and kernel parameter for data set 2
[Ks,Kt] = conKernelMatrix(train_data',test_data',kernel1,kernel2,lrank);

temp = train_label;
temp(temp==0) = 2;
def_n = length(find(temp==1));
nodef_n = length(find(temp==2));
c = 2; % classes
cost(1) = 1; 
cost(2) = nodef_n/def_n;
Css = 0;
Cst = 0;
for i=1:c
    idx = temp == i;
    ksi = Ks(idx,:);
    Css = Css+cost(i)*(ksi'*ksi)+reg*eye(size(ksi,2));
    
    dist = pdist2(ksi,Kt);
    dist = exp(-dist);
    Cst = Cst+cost(i)*ksi'*dist*Kt;
end

Ctt = Kt'*Kt+reg*eye(size(Kt,2));

[Ws,Wt,s] = eigDecomposition(Css,Ctt,Cst);

dim = ceil(size(train_data,1)*0.15); % tune the projected dimension 

Wxx = Ws(:,1:dim);
Wyy = Wt(:,1:dim);
train_new = Ks*Wxx;
test_new = Kt*Wyy;

%% LR 
model = train(train_label', sparse(real(train_new)),'-s 0 -c 1 -B -1 -q'); % num * fec
[predict_label,~, prob_estimates] = predict(test_label', sparse(real(test_new)), model,'-b 1');
score = prob_estimates(:,1)';
predict_label = predict_label';
