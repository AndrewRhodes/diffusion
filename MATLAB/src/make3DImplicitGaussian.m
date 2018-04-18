% % Andrew Rhodes
% December 2017
% make3DGaussianMatrix
%
% Builds the 3D discrete Gaussian
%
% input[]:
% input[]:
% input[]:
% input[]:
% input[]:
% input[]:
% input[]:
%
% output[G]:



function G = make3DImplicitGaussian(x, y, z, sigma, spacing, band, numsigmas, LimitFarPoints)

% Prob3D = chi2cdf(numsigma^2,dimension)

% Preferable that sigma < spacing

% sigma = 0.25;
% spacing = 0.5;


% if sigma == spacing
%     SpacingEqualSigma = true;
% end


% clc
% numsigmas = 4;
% sigma = 0.1
% spacing = 0.3
% LimitFarPoints=1

TotalSigma = numsigmas * sigma;
SpacingSigmaRatio =  ceil(TotalSigma / spacing);


vec1d = -SpacingSigmaRatio*spacing:spacing:SpacingSigmaRatio*spacing;


[xg,yg,zg] = meshgrid(vec1d, vec1d, vec1d);

weights = exp( - (xg.^2 + yg.^2 + zg.^2) ./ (2*sigma^2) );
weights = reshape(weights, [],1);
weights = weights ./ sum(weights(:));


PTS = [reshape(xg,[],1),reshape(yg,[],1),reshape(zg,[],1)]./spacing;
PTS = double(int16(PTS));

if LimitFarPoints
    TooFarPts = reshape(sqrt(xg.^2 + yg.^2 + zg.^2) > TotalSigma,[],1);
    
    weights(TooFarPts) = [];
    PTS(TooFarPts,:) = [];
end


% Possibly could reduce the size of weights to only include points within a
% ball of radius 4*spacing, if sigma<spacing


Nx = length(x);
Ny = length(y);
Nz = length(z);


StencilSize = length(weights);

Gi = repmat((1:length(band))',1, StencilSize);
Gj = zeros(size(Gi));
Gs = zeros(size(Gi));


[j,i,k] = ind2sub([Ny,Nx,Nz], band);

tic
for c = 1 : StencilSize
    
    ii = i + PTS(c,2);
    jj = j + PTS(c,1);
    kk = k + PTS(c,3);
    Gs(:,c) = weights(c);
    
    Gj(:,c) = sub2ind([Ny,Nx,Nz],jj,ii,kk);
%     Gj(:,c) = sub2ind([Nx,Ny,Nz],ii,jj,kk);

    
end
toc


G = sparse(Gi(:), Gj(:), Gs(:), length(band), Nx*Ny*Nz, nnz(Gs));


% Gout = G(:, setdiff(1:(Nx*Ny*Nz),Band));

% nnzG = nnz(G);



% if nnzG ~= nnz(G)
%     disp('Lost some non-zero coefficients (from outside the outerband)')
% end


% Normalize the rows to unity after simplifying by band
G = G(:,band);

Glength = length(G);
% 
G = sparse(1:Glength,1:Glength, 1./sum(G,2)) * G;



end













