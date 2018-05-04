







function [ItL, Eplot, CP, CPFACE] = makeImplicitLaplaceBeltrami( PointCloudIn, porder, Lorder, spacing, dim, OrderBDF, StepSize, Alpha)

if Alpha < 0
    error('Expected Alpha >= 0, but input is ''%0.3f''', Alpha)
end

if ~isfield(PointCloudIn, 'Location')
    error('PointCloudIn must contain the field ''Location''.')
end

if ~isfield(PointCloudIn, 'Face')
    error('PointCloudIn must contain the field ''Face''.')
end

if ~isfield(PointCloudIn, 'Signal')
    error('PointCloudIn must contain the field ''Signal''.')
end



bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));


MinPoint = min(PointCloudIn.Location) - bandwidth - 2*spacing;
MaxPoint = max(PointCloudIn.Location) + bandwidth + 2*spacing;


x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';



% [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloudIn.Face, PointCloudIn.Location, spacing, MinPoint, porder, Lorder/2);
[IJK,DD,CP,XYZ,CPFACE] = tri2cp_helper(spacing, MinPoint, bandwidth, PointCloudIn.Face, PointCloudIn.Location, 1, 10);
% Note: DD is squared distance
DIST = sqrt(DD);


% Notice the order of j,i,k, and y,x,z
BandSearchSize = [length(y1d), length(x1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,2), IJK(:,1), IJK(:,3));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matric Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Lc = laplacian_3d_matrix(x1d, y1d, z1d, Lorder, Band);

Eplot = interp3_matrix(x1d, y1d, z1d, PointCloudIn.Location(:,1), PointCloudIn.Location(:,2), PointCloudIn.Location(:,3), porder, Band);
% Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);

Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);
% Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);

L = lapsharp(Lc, Ecp);


[Lheight, Lwidth] = size(L);




ItL = cell(OrderBDF,1);

ItL{1,1} = speye(Lheight, Lwidth) - Alpha * StepSize * L;


if OrderBDF == 2 || OrderBDF == 3 || OrderBDF == 4
%     I23tL
    ItL{2,1} = speye(Lheight, Lwidth) - (2/3) * Alpha * StepSize * L;
end

if OrderBDF == 3 || OrderBDF == 4
%     I611tL
    ItL{3,1} = speye(Lheight, Lwidth) - (6/11) * Alpha * StepSize * L;
end

if OrderBDF == 4
%     I1225tL
    ItL{4,1} = speye(Lheight, Lwidth) - (12/25) * Alpha * StepSize * L;
end









end