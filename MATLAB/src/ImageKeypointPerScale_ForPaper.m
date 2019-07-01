% Andrew Rhodes
% ASEL
% December 2018


close all
clear
clc

global ProjectRoot;

UseNMS = 0;
alpha = 1;

FileLocation = strcat(ProjectRoot,'/models/image/');
PointCloudLocation = '/home/andrew/GitAndrew/diffusion/MATLAB/models/image/';
PointCloudName = {'Voorhees_e1_90688', 'Woodburn_e1_141198'};
KeypointLocation = strcat(ProjectRoot, '/data/imagekeypointdata/')
 
 
TitleNames = {'SIFT', 'Mesh Euc.', 'Cotangent', 'Umbrella'};


modelFolder = {'Voorhees','Woodburn'}; %  {'Voorhees','Woodburn','Woodburn2'}; % 

LBO = {'SIFT', 'mesh', 'cot', 'umb'};
 
LBOLength = length(LBO);
modelFolderLength = length(modelFolder);


k = sqrt(2);
l_scale = 1/k;
set_ltau = @(e_bar) 2*e_bar;
set_t0 = @(e_bar) e_bar ;

LineStyles = {'b-', 'r--', 'g-.', 'k-', 'c-'};


for jj = 1 : LBOLength
    
    
    
    figure
    
           
    maxx = zeros(modelFolderLength, 1);
    
    for j = 1 : modelFolderLength               
        
        
        if strcmpi(modelFolder{j}, 'Voorhees')

            ImageName = 'VoorheesComputingCenter_831x436.jpg';
            Image = double(rgb2gray(imread(strcat(FileLocation,ImageName))));
            Image = imresize(Image,0.5);
            Image = Image./255;
            
            N = 9;
        elseif strcmpi(modelFolder{j}, 'Woodburn')

            ImageName = 'WVUWoodburn_524x358.jpg';
            Image = double(rgb2gray(imread(strcat(FileLocation,ImageName))));
            Image = Image./255;

            N = 10;            
        elseif strcmpi(modelFolder{j}, 'Woodburn2')

            ImageName = 'WVUWoodburn_438x284.jpg';
            Image = double(rgb2gray(imread(strcat(FileLocation,ImageName))));
            Image = Image./255;    

            N = 10;            
        end
        
        
        ImageSize = size(Image);
        NumVertex = numel(Image);
        NumPixels = numel(Image);

        
        [PointCloud.Location,PointCloud.Face] = read_off( strcat(PointCloudLocation, PointCloudName{j}, '.off') );
        load( strcat(PointCloudLocation, PointCloudName{j}, '_Neighbors.mat'),  'Neighbors' )
        PointCloud.Location = PointCloud.Location';
        PointCloud.Face = PointCloud.Face';        
        PointCloud.LocationCount = length(PointCloud.Location);
        PointCloud.FaceCount = length(PointCloud.Face);
        PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
        PointCloud.Signal = reshape(Image, [], 1);
        PointCloud = findMeshResolution(PointCloud, 'Model');

        
        l_tau = set_ltau(PointCloud.Resolution);
        t_0 = set_t0(PointCloud.Resolution);
        
        
        load( strcat(KeypointLocation, modelFolder{j}, '/Keypoint_', LBO{jj}, '.mat') , strcat('Keypoint_', LBO{jj}) )
        
        if jj == 1 % SIFT
            Keypoint = Keypoint_SIFT;
        elseif jj == 2 % mesh
            Keypoint = Keypoint_mesh;
        elseif jj == 3 % cot
            Keypoint = Keypoint_cot;
        elseif jj == 4 % umb
            Keypoint = Keypoint_umb;
        end
 
        
        
        
        
        if jj == 1 % SIFT
            
            [HistHeight, HistEdge, HistBin] = histcounts(Keypoint(3,:), N, 'BinWidth', 1);

%             plot(HistEdge(1:end-1), HistHeight/max(HistHeight), LineStyles{j}, 'Linewidth', 4)
            semilogy(HistEdge(1:end-1), HistHeight, LineStyles{j}, 'Linewidth', 4)
            hold on
            maxx(j) = max(Keypoint(3,:));
                
        else
            
            LevelTable = tabulate(Keypoint.Level);
            
            [tn, tau_n, sigma_n] = findScaleStep(k, t_0, alpha, N);                    
            
            if jj == 2 % mesh
                
%                 plot(sigma_n(LevelTable(:,1)), ...
%                 LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{j}, 'Linewidth', 4)
                semilogy(sigma_n(LevelTable(:,1)), ...
                    LevelTable(:,2), LineStyles{j}, 'Linewidth', 4)
                hold on
                maxx(j) = max(sigma_n(LevelTable(:,1))./PointCloud.Resolution);
            
            elseif jj == 3 % cot
                
%                 plot(sigma_n(LevelTable(:,1)), ...
%                 LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{j}, 'Linewidth', 4)
            
                semilogy(sigma_n(LevelTable(:,1)), ...
                    LevelTable(:,2), LineStyles{j}, 'Linewidth', 4)
                hold on
                maxx(j) = max(sigma_n(LevelTable(:,1))./sqrt(PointCloud.Resolution));
            
            elseif jj == 4 % umb
                
%                 plot(sigma_n(LevelTable(:,1)), ...
%                 LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{j}, 'Linewidth', 4)
                
                   semilogy(sigma_n(LevelTable(:,1)), ...
                    LevelTable(:,2), LineStyles{j}, 'Linewidth', 4)
                hold on
                maxx(j) = max(PointCloud.Resolution*sigma_n(LevelTable(:,1)));
                
            end
               
        end 
                        
        drawnow
        hold on
        
        clear Keypoint PointCloud LevelTable tau
    
    end
        
    
    title( TitleNames{jj})
%     hold on
    xlabel({'$\sigma$'}, 'interpreter','Latex')
%     if (jj == 1) || (jj == 4)
%         xlabel({'$\sigma$'}, 'interpreter','Latex')
%     elseif jj == 3
%         xlabel({'$\sigma / \sqrt{\bar{e}}$'}, 'interpreter','Latex')
%     else
%         xlabel({'$\sigma / \bar{e}$'}, 'interpreter','Latex')
%     end
    ylabel('No. Keypoints')
    ax = gca;
    ax.XAxis.FontSize = 50;
    ax.XLabel.FontSize = 75;
    ax.YAxis.FontSize = 50;
    ax.Title.FontSize = 50;
    
    
%     xlim([1 10])
    xlim([1 min(maxx)])
    xlim([2 15])
    title([])
    
%     sprintf('%s: mean: %0.4f   std: %0.4f', TitleNames{j}, mean(x), std(x))
     
%     pause
end





















