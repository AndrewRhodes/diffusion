% Andrew Rhodes
% West Virginia University
% 
% Computer the symmetric Gaussian matrix from mesh
%
% input[filename]: '.off' filename of triangle mesh
% input[opt]: program options.
%             opt.dtype: distance calculation method.
%             dtype = 'euclidean' or 'geodesic';
%             Default: 'euclidean'
%           
%             opt.atype: use area weighting
%             atype = false or true
%             Default: false
%           
%             opt.rho: cutoff for Gaussian function evaluation
%             Default: 3, must > 0
% 
%             opt.hs: Scaling factor for neighborhood size to the parameter
%                     h, where h^2 = 4t
%             Default: 2, must > 0
%             
%
% output[G]: Gaussian weighting matrix
% output[h]: Gaussian width: h^2 = 4t


% Main difference from symmshlp_matrix, is that we do NOT consider triangle
% areas. To Do: Add triangle areas into another code.



% function [G, h] = mesh_Gauss(filename, opt)
function W = mesh_Gauss(filename, opt)


if nargin < 1
    error('Too few input arguments');	 
end

opt = parse_opt(opt);

if opt.hs <= 0 || opt.rho <= 0
	error('Invalid values in opt');
end

if opt.atype
    [I, J, S, A, h] = meshGauss(filename, opt);
else
    [I, J, S, h] = meshGauss(filename, opt);
    
    W = sparse(I,J,S);
    Wsum = sparse(1:length(W), 1:length(W), 1./sum(W,2));    
    W = Wsum * W;

end



% G = sparse(I, J, S);


end



function opt = parse_opt(opt)

if ~isfield(opt, 'hs')
	opt.hs = 2;
end


if ~isfield(opt, 'rho')
	opt.rho = 3;
end

if ~isfield(opt, 'htype')
% 	opt.htype = 'ddr';
    opt.htype = 0;
else
    if strcmp(opt.htype, 'ddr')
        opt.htype = 0;
    elseif strcmp(opt.htype, 'ppr')
        opt.htype = 1;  
    end
end

if ~isfield(opt, 'dtype')
% 	opt.dtype = 'euclidean';
    opt.dtype = 0;
else 
    if strcmp(opt.dtype, 'euclidean')
        opt.dtype = 0;
    elseif strcmp(opt.dtype, 'geodesic')
        opt.dtype = 1;  
    end
end

if ~isfield(opt, 'atype')
% 	opt.atype = 0;
    opt.atype = false;
end

end




