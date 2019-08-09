function K = km_kernel(X1,X2,ktype,kpar)
% KM_KERNEL calculates the kernel matrix between two data sets.
% Input:	- X1, X2: data matrices in row format (data as rows)
%			- ktype: string representing kernel type
%			- kpar: vector containing the kernel parameters
% Output:	- K: kernel matrix
% USAGE: K = km_kernel(X1,X2,ktype,kpar)
%
% Author: Steven Van Vaerenbergh (steven *at* gtas.dicom.unican.es), 2010.
% Id: km_kernel.m v1.0
% This file is part of the Kernel Methods Toolbox (KMBOX) for MATLAB.
% http://sourceforge.net/p/kmbox
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, version 3 (as included and available at
% http://www.gnu.org/licenses).

switch ktype
	case 'gauss'	% Gaussian kernel
		% this calculation introduces a small numerical error w.r.t. to the
		% direct calculation, but it is much faster.
		sgm = kpar;	% kernel width
		
		dim1 = size(X1,1);
		dim2 = size(X2,1);
		
% 		if (dim1*dim2 >= 1E6)	% warning for large computations
% 			fprintf('Warning: computation of %d kernel elements might be slow. Press key to continue?\n',dim1*dim2);
% 			pause
% 		end
		
		norms1 = sum(X1.^2,2);
		norms2 = sum(X2.^2,2);
		
		mat1 = repmat(norms1,1,dim2);
		mat2 = repmat(norms2',dim1,1);
		
		distmat = mat1 + mat2 - 2*X1*X2';	% full distance matrix
		K = exp(-distmat/(2*sgm^2));
		
	case 'gauss-diag'	% only diagonal of Gaussian kernel
		sgm = kpar;	% kernel width
		K = exp(-sum((X1-X2).^2,2)/(2*sgm^2));
	
    case 'linear'
        K = X1*X2';
        
    case 'linear-diag'
        K = X1*X2';
        K = diag(K);
        
    case 'poly'	% polynomial kernel
		p = kpar(1);	% polynome order
		c = kpar(2);	% additive constant
		
		K = (X1*X2' + c).^p; 
        
	case 'poly-diag'	% polynomial kernel
		p = kpar(1);	% polynome order
		c = kpar(2);	% additive constant
		
		K = (X1*X2' + c).^p;
        K = diag(K);
		
	otherwise	% default case
		error ('unknown kernel type')
end