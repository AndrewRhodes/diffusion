% Andrew Rhodes
% ASEL
% December 2018


close all
clear
clc

global ProjectRoot;

UseNMS = 1;
t_scale = 1/sqrt(2);
level_min = 3;
t_range = 1/2;


LBO = {'Mesh', 'Mesh_euc', 'Umbrella', 'Cotangent'};
modelFolder = {'bunny','armadillo','buddha','dragon','itokawa'}; % 
modelNames = {'Bunny_e1','Armadillo_e1','Buddha_e1', 'Dragon_e1', 'Itokawa_e1'}; % 

LBOLength = length(LBO);
modelFolderLength = length(modelFolder);
modelNamesLength = length(modelNames);


SamplePercVecBunny = [1, 0.98, 0.95, 0.93, 0.9, 0.88, 0.85, 0.83]; % bunny
SamplePercVecArmadillo = [0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3]; % armadillo
SamplePercVecBuddha = [0.065, 0.06, 0.058, 0.055, 0.053, 0.05, 0.048, 0.045]; % buddha
SamplePercVecDragon = [0.072, 0.07, 0.068, 0.065, 0.063, 0.06, 0.058, 0.055]; % dragon
SamplePercVecItokawa = [0.65, 0.62, 0.6, 0.58, 0.55, 0.52, 0.5, 0.48]; % Itokawa


SamplePercVec = [SamplePercVecBunny; SamplePercVecArmadillo; ...
                SamplePercVecBuddha; SamplePercVecDragon; ...
                SamplePercVecItokawa];

FileLocationWD = '/media/andrew/WDRhodes/diffusiondata/';

Error = cell(modelNamesLength, 7, LBOLength);
NoMatch = cell(modelNamesLength, 7, LBOLength);


for kkk = 1 : LBOLength
    
    for k = 1 : modelNamesLength
        %     sprintf('Model: %s, LBO: %s', modelFolder{k}, LBO{kkk})
        sprintf('Model: %s', modelFolder{k})
        
        
        
        ModelFolder = strcat(modelFolder{k},'/');
        FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');
        ListOriginal = dir(strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Sample/',LBO{kkk},'/'));
        ListCases = dir(strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Sample/',LBO{kkk},'/Cases/'));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pick the correct 'Truth' data
        for i = 1 : length([ListOriginal.isdir])
            if UseNMS
                IsThisCorrect = regexp(ListOriginal(i).name, 'NMS', 'start');
            else
                IsThisCorrect = regexp(ListOriginal(i).name, 'Keypoint', 'start');
            end
            if IsThisCorrect == 1
                TruthIndex = i;
                SplitName = strsplit(ListOriginal(i).name,'_');
                FaceSizesStrOrig = SplitName{3};
                SplitStepsOriginal = strsplit(SplitName{4},{'N','.'});
                NumStepsOriginal = str2double(SplitStepsOriginal{2});
                break;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the point cloud original
        FaceSizesOrig = str2num(FaceSizesStrOrig);
        ModelOriginal = strcat(modelNames{k},'_',FaceSizesStrOrig);
        [PointCloudOriginal.Location, PointCloudOriginal.Face, PointCloudOriginal.Normal, PointCloudOriginal.Signal]...
            = read_ply_all_elements( strcat(ProjectRoot,'/models/object/',ModelFolder,'Sample/',ModelOriginal,'.ply') );
        
        PointCloudOriginal.LocationCount = size(PointCloudOriginal.Location,1);
        PointCloudOriginal.FaceCount = size(PointCloudOriginal.Face, 1);
        PointCloudOriginal.FaceArea = findFaceArea(PointCloudOriginal.Location,PointCloudOriginal.Face);
        PointCloudOriginal = findMeshResolution(PointCloudOriginal, 'Model');
        
        FileNameNeighbors = strcat(ModelOriginal,'_Neighbors.mat');
        if ~exist( strcat( FileLocationNeighbors, FileNameNeighbors), 'file')
            [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
            save(strcat( FileLocationNeighbors, FileNameNeighbors) ,'Neighbors', '-v7.3')
        else
            load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
        end
        PointCloudOriginal = findLocalResolution(PointCloudOriginal, Neighbors.Connect);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load Keypoint Original
        if UseNMS
            load(strcat(ListOriginal(TruthIndex).folder, '/', ListOriginal(TruthIndex).name));
            KeypointOriginal = NMSKeypoint;
            clear NMSKeypoint
        else
            load(strcat(ListOriginal(TruthIndex).folder, '/', ListOriginal(TruthIndex).name));
            KeypointOriginal = Keypoint;
            clear Keypoint
        end
        
        
        
        % Place different LBO types here
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Select and order the cases
        FaceSizesNum = [];
        CorrectFiles = [];
        NumStepsStr = [];
        Indices = [];
        for i = 1 : length([ListCases.isdir])
            if UseNMS
                IsThisCorrect = regexp(ListCases(i).name, 'NMS', 'start');
            else
                IsThisCorrect = regexp(ListCases(i).name, 'Keypoint', 'start');
            end
            if IsThisCorrect == 1
                SplitName = strsplit(ListCases(i).name,'_');
                FaceSizesNum = [FaceSizesNum; str2num(SplitName{3})];
                CorrectFiles = [CorrectFiles; i];
                
                SplitSteps = strsplit(SplitName{4},{'N','.'});
                NumStepsStr = [NumStepsStr; SplitSteps{2}];
                
                Indices = [Indices; i];
            end
        end
        NumSteps = str2num(NumStepsStr);
        %     FaceSizesNum = str2num(FaceSizesStr);
        [SortVal, SortOrd] = sort(FaceSizesNum, 'descend');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Repeatibility between original and cases
        for kk = 1 : length(CorrectFiles)
            
            clear PointCloud
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load the point cloud case
            ModelIter = strcat(modelNames{k},'_',num2str(FaceSizesNum(SortOrd(kk))));
            [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
                = read_ply_all_elements( strcat(ProjectRoot,'/models/object/',ModelFolder,'Sample/',ModelIter,'.ply') );
            
            PointCloud.LocationCount = size(PointCloud.Location,1);
            PointCloud.FaceCount = size(PointCloud.Face, 1);
            PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
            PointCloud = findMeshResolution(PointCloud, 'Model');
            
            FileNameNeighbors = strcat(ModelIter,'_Neighbors.mat');
            if ~exist( strcat( FileLocationNeighbors, FileNameNeighbors), 'file')
                [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
                save(strcat( FileLocationNeighbors, FileNameNeighbors) ,'Neighbors', '-v7.3')
            else
                load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
            end
            PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load Keypoint Case
            if UseNMS
                load(strcat(ListCases(1).folder, '/', ListCases(Indices(SortOrd(kk))).name))
                KeypointIter = NMSKeypoint;
                clear NMSKeypoint
            else
                load(strcat(ListCases(1).folder, '/', ListCases(Indices(SortOrd(kk))).name))
                KeypointIter = Keypoint;
                clear Keypoint
            end
            
            
            [Error{k, kk, kkk}, NoMatch{k, kk, kkk}] = findKeypointErrorSample(PointCloudOriginal, KeypointOriginal,...
                NumStepsOriginal, PointCloud, KeypointIter,...
                NumSteps(SortOrd(kk)), level_min, t_scale, t_range);
            
        end
        
        
        
        
    end
    
end


positions = [1, 1.1429, 1.2857, 1.4286, 1.5714, 1.7143, 1.8571];
% LineStyles = {'k-', 'r--', 'b:'};
LineStyles = {'k-', 'g--', 'r:', 'b-.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for jjj = 1 : LBOLength
    for jj = 1 : modelNamesLength
        cp = positions + (jj - 1)*2;
        Data = zeros(size(Error,2),1);
        for j = 1 : size(Error,2)
            Data(j) = Error{jj, j, jjj}.ScaleRepeat;
        end        
        p(jjj) = plot(cp, Data, LineStyles{jjj}, 'Linewidth', 4);
    end
end

axis([0,2*modelNamesLength,0,1])
ylim([0 1])
%     xlabel('Model')
ylabel('Percentage')
ax = gca;
ax.YTick = [0,0.2,0.4,0.6,0.8,1];
ax.XTick = 1.25 + ((1:modelNamesLength)-1)*2;
ax.XTickLabel = modelFolder;
xtickangle(ax,-30)
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;
% legend(p, {'Mesh Geo.','Umb.','Cot.'}, 'FontSize', 30)
legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal', 'Location', 'northoutside')
% legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for jjj = 1 : LBOLength
    for jj = 1 : modelNamesLength
        cp = positions + (jj - 1)*2;
        Data = zeros(size(Error,2),1);
        for j = 1 : size(Error,2)
            Data(j) = Error{jj, j, jjj}.RelativeRepeat;
        end        
        p(jjj) = plot(cp, Data, LineStyles{jjj}, 'Linewidth', 4);
    end
end

axis([0,2*modelNamesLength,0,1])
ylim([0 1])
%     xlabel('Model')
ylabel('Percentage')
ax = gca;
ax.YTick = [0,0.2,0.4,0.6,0.8,1];
ax.XTick = 1.25 + ((1:modelNamesLength)-1)*2;
ax.XTickLabel = modelFolder;
xtickangle(ax,-30)
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;
% legend(p, {'Mesh Geo.','Umb.','Cot.'}, 'FontSize', 30)
legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal', 'Location', 'northoutside')
% legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for jjj = 1 : LBOLength
    for jj = 1 : modelNamesLength
        cp = positions + (jj - 1)*2;
        Data = zeros(size(Error,2),1);
        for j = 1 : size(Error,2)
            Data(j) = Error{jj, j, jjj}.Count;
        end        
        p(jjj) = plot(cp, Data, LineStyles{jjj}, 'Linewidth', 4);
        ylims(jj,:) = get(gca, 'YLim');
    end
end


xlim([0,2*modelNamesLength])
set(gca, 'Ylim', [0 max(ylims(:,2))]);
ylabel('Count')
ax = gca;
%         ax.YTick = [0,40,80,120,160];
ax.XTick = 1.25 + ((1:modelNamesLength)-1)*2;
ax.XTickLabel = modelFolder;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;
xtickangle(ax,-30)
% legend(p, {'Mesh Geo.','Umb.','Cot.'}, 'FontSize', 30)
legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal', 'Location', 'northoutside')
% legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal')
















