function [Ws,Wt,s] = eigDecomposition(Css,Ctt,Cst)
% Here, we use the singular value decompose (SVD) as paper for eigenvalue decomposition.
% 
% Or you can use generalized eigenvalue decompostion to solve it. 
% [V,D] = eig(A,B) returns diagonal matrix D of generalized eigenvalues and full matrix V 
% whose columns are the corresponding right eigenvectors, so that A*V = B*V*D
% 
% In this paper, A = [0, Cst; Cst', 0]; B = [Css 0; 0 Ctt]. 
%

[u,s,v] = svd(Css);
X = u*s^(-0.5)*v';
[u,s,v] = svd(Ctt);
Y = u*s^(-0.5)*v';

H = X*Cst*Y;
[u,s,v] = svd(H);
Ws = X*u;
Wt = Y*v;

Ws = Ws./repmat(sqrt(sum(abs(Ws).^2)),size(Ws,1),1);
Wt = Wt./repmat(sqrt(sum(abs(Wt).^2)),size(Wt,1),1);

s = diag(s);


