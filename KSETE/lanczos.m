function [ T, V ] = lanczos( A, m)
% lanczos Summary of this function goes here
%   Detailed explanation goes here 
% INPUTS:
%   (1) A - a symmetric matrix of size n*n.
%   (2) m - optionally a number of iterations m (m<=n)(as default, let m=n).
% OUTPUTS:
%   (1) T - a tridiagonal real symmetric matrix T=V'AV of size m*m.
%   (2) V - an n*m matrix V with orthonormal columns.


if nargin==1
    m = size(A,1);
end

% An arbitrary vector with euclidean norm 1.
v = zeros(size(A,1),1);
ind = randperm(size(A,1),size(A,1));
v(ind(1)) = 1;
ind_count = 1;

V = zeros(size(A,1), m);

% 
alpha = zeros(m,1);
beta = zeros(m,1);

for i=1:m
    
    if i==1
        w = A*v; % w is a 
        alpha(i) = w'*v; % scalar product, dot product
        w = w - alpha(i)*v;
        V(:,i) = v;
    else
        beta(i) = norm(w);
        % Update v
        if beta(i)~=0
            v = w/beta(i);
            V(:,i) = v;
        else
            v = zeros(size(A,1),1);
            ind_count = ind_count + 1;
            v(ind(ind_count)) = 1; % an arbitrary vector with Euclidean norm 1 and it is orthogonal to all of V(:,1),...,V(:,i-1).
            V(:,i) = v;
        end
        
        w = A*v; % 
        alpha(i) = w'*v; % scalar product, dot product
        w = w - alpha(i)*v - beta(i)*V(:,i-1);        
    end
end

T = V'*A*V;

end

