










function [KeypointA, KeypointB, KeypointC, KeypointD] = findAllKeypoints(PointCloud, DoG, ScaleParameter, Neighbors)

NumDoGLevels = size(DoG,2);

KeypointA.Location = zeros(PointCloud.LocationCount,1);
KeypointA.Level = zeros(PointCloud.LocationCount,1);
KeypointA.Scale = zeros(PointCloud.LocationCount,1);
KeypointA.Sign = zeros(PointCloud.LocationCount,1);

KeypointB.Location = zeros(PointCloud.LocationCount,1);
KeypointB.Level = zeros(PointCloud.LocationCount,1);
KeypointB.Scale = zeros(PointCloud.LocationCount,1);
KeypointB.Sign = zeros(PointCloud.LocationCount,1);

KeypointC.Location = zeros(PointCloud.LocationCount,1);
KeypointC.Level = zeros(PointCloud.LocationCount,1);
KeypointC.Scale = zeros(PointCloud.LocationCount,1);
KeypointC.Sign = zeros(PointCloud.LocationCount,1);

KeypointD.Location = zeros(PointCloud.LocationCount,1);
KeypointD.Level = zeros(PointCloud.LocationCount,1);
KeypointD.Scale = zeros(PointCloud.LocationCount,1);
KeypointD.Sign = zeros(PointCloud.LocationCount,1);


% KeypointE.Location = 0;
% KeypointE.Level = 0;
% KeypointE.Scale = 0;
% KeypointE.Sign = 0;
% KeypointE.DoG = 0;

isKeypointMethodA = 0;
isKeypointMethodB = 0;
isKeypointMethodC = 0;
isKeypointMethodD = 0;
% isKeypointMethodE = 0;
NumKeypointA = 0;
NumKeypointB = 0;
NumKeypointC = 0;
NumKeypointD = 0;
% NumKeypointE = 0;
isSignA = 0;
isSignB = 0;
isSignC = 0;
isSignD = 0;
% isSignE = 0;

WaitBar = waitbar(0, sprintf('Checking Vertex %i of %i', 0, PointCloud.LocationCount));

for i = 1 : PointCloud.LocationCount
    
    
    
%     SignDoG = movprod(sign(DoG(i,2:end)),2)
%     SignDiffDoG = movprod(sign(diff(DoG(i,2:end))),2)
%     
%     SlopeChangesDoG = find(SignDiffDoG < 0)
%     SlopeChangesDoG(SlopeChangesDoG==1) = []
%     SlopeChangesDoG(SlopeChangesDoG==length(SignDiffDoG)) = []
%     
% 
%     for k = 1 : length(SlopeChangesDoG)
%         Neighbors = DoG(i,[SlopeChangesDoG(k)-1,SlopeChangesDoG(k)+2]);
%         (Neighbors(1) -2* DoG(i,SlopeChangesDoG(k)+1) + Neighbors(2))
%         
%         if DoG(i,SlopeChangesDoG(k)+1) < 0
%             if all(DoG(i,SlopeChangesDoG(k)+1) < Neighbors)
%                 isKeypointMethodE = 1;
%                 isSignE = -1;
%             end
%         elseif DoG(i,SlopeChangesDoG(k)+1) > 0
%             if all(DoG(i,SlopeChangesDoG(k)+1) > Neighbors)
%                 isKeypointMethodE = 1;
%                 isSignE = 1;
%             end
%         end
%         
%         if isKeypointMethodE
%             NumKeypointE = NumKeypointE + 1;
%             KeypointE.Location(NumKeypointE, 1) = i;
%             KeypointE.Level(NumKeypointE, 1) = SlopeChangesDoG(k)+1;
%             KeypointE.Scale(NumKeypointE, 1) = ScaleParameter(SlopeChangesDoG(k)+1) + PointCloud.ResolutionLocal(i);
%             KeypointE.Sign(NumKeypointE, 1) = isSignE;
%             KeypointE.DoG(NumKeypointE, 1) = DoG(i,SlopeChangesDoG(k)+1);
%             
%             isSignE = 0;
%             isKeypointMethodE = 0;
%         end
%         
%     end
    
        
    CurrentNeighborsConnect = Neighbors.Connect{i,1};
    CurrentNeighborsDistance = Neighbors.Distance{i,1};
    
    
    for j = 2 : NumDoGLevels - 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract DoG Values for consideration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        DoGCurrent = DoG(i,j);
        DoGAbove = DoG(i, j+1);
        DoGBelow = DoG(i, j-1);
        
        NeighborValuesBelowConnect = DoG(CurrentNeighborsConnect, j-1);
        NeighborValuesCurrentConnect = DoG(CurrentNeighborsConnect, j);
        NeighborValuesAboveConnect = DoG(CurrentNeighborsConnect, j+1);
        
        NeighborValuesBelowDistance = DoG(CurrentNeighborsDistance, j-1);
        NeighborValuesCurrentDistance = DoG(CurrentNeighborsDistance, j);
        NeighborValuesAboveDistance = DoG(CurrentNeighborsDistance, j+1);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compare the DoG Values to the neighbors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Connectivity Neighbors
        if DoGCurrent > DoGAbove
            if DoGCurrent > DoGBelow
                if all(DoGCurrent > NeighborValuesCurrentConnect)
                    isKeypointMethodC = 1;
                    isSignC = 1;
                    if all(DoGCurrent > NeighborValuesBelowConnect)
                        if all(DoGCurrent > NeighborValuesAboveConnect)
                            isKeypointMethodA = 1;
                            isSignA = 1;
                        end
                    end
                end
            end
        elseif DoGCurrent < DoGAbove
            if DoGCurrent < DoGBelow
                if all(DoGCurrent < NeighborValuesCurrentConnect)
                    isKeypointMethodC = 1;
                    isSignC = -1;
                    if all(DoGCurrent < NeighborValuesBelowConnect)
                        if all(DoGCurrent < NeighborValuesAboveConnect)
                            isKeypointMethodA = 1;
                            isSignA = -1;
                        end
                    end
                end
            end
        end
        
                    
                    
        % Distance Neighbors
        if DoGCurrent > DoGAbove
            if DoGCurrent > DoGBelow
                if all(DoGCurrent > NeighborValuesCurrentDistance)
                    isKeypointMethodD = 1;
                    isSignB = 1;
                    if all(DoGCurrent > NeighborValuesBelowDistance)
                        if all(DoGCurrent > NeighborValuesAboveDistance)
                            isKeypointMethodB = 1;
                            isSignD = 1; 
                        end
                    end
                end
            end
        elseif DoGCurrent < DoGAbove
            if DoGCurrent < DoGBelow
                if all(DoGCurrent < NeighborValuesCurrentDistance)
                    isKeypointMethodD = 1;
                    isSignB = -1;
                    if all(DoGCurrent < NeighborValuesBelowDistance)
                        if all(DoGCurrent < NeighborValuesAboveDistance)
                            isKeypointMethodB = 1;
                            isSignD = -1;
                        end
                    end
                end
            end
        end
        
        
        
        if isKeypointMethodA
            NumKeypointA = NumKeypointA + 1;
            KeypointA.Location(NumKeypointA, 1) = i;
            KeypointA.Level(NumKeypointA, 1) = j;
            KeypointA.Scale(NumKeypointA, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i);
            KeypointA.Sign(NumKeypointA, 1) = isSignA;
            
            isSignA = 0;
            isKeypointMethodA = 0;
        end
        
        
        
        if isKeypointMethodB
            NumKeypointB = NumKeypointB + 1;
            KeypointB.Location(NumKeypointB, 1) = i;
            KeypointB.Level(NumKeypointB, 1) = j;
            KeypointB.Scale(NumKeypointB, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i);
            KeypointB.Sign(NumKeypointB, 1) = isSignB;
            
            isSignB = 0;
            isKeypointMethodB = 0;
        end
        
        
        
        if isKeypointMethodC
            NumKeypointC = NumKeypointC + 1;
            KeypointC.Location(NumKeypointC, 1) = i;
            KeypointC.Level(NumKeypointC, 1) = j;
            KeypointC.Scale(NumKeypointC, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i);
            KeypointC.Sign(NumKeypointC, 1) = isSignC;
            
            isSignC = 0;
            isKeypointMethodC = 0;
        end
        
        
        
        if isKeypointMethodD
            NumKeypointD = NumKeypointD + 1;
            KeypointD.Location(NumKeypointD, 1) = i;
            KeypointD.Level(NumKeypointD, 1) = j;
            KeypointD.Scale(NumKeypointD, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i);
            KeypointD.Sign(NumKeypointD, 1) = isSignD;
            
            isSignD = 0;
            isKeypointMethodD = 0;
        end
        
        
    end
    
    
    waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Checking Vertex %i of %i', i, PointCloud.LocationCount));
    
end




waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Checking Vertex for Keypoint Complete'));
close(WaitBar)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check every vertex across scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZeroLogicA = KeypointA.Location==0;
KeypointA.Location(ZeroLogicA) = [];
KeypointA.Level(ZeroLogicA) = [];
KeypointA.Scale(ZeroLogicA) = [];
KeypointA.Sign(ZeroLogicA) = [];

ZeroLogicB = KeypointB.Location==0;
KeypointB.Location(ZeroLogicB) = [];
KeypointB.Level(ZeroLogicB) = [];
KeypointB.Scale(ZeroLogicB) = [];
KeypointB.Sign(ZeroLogicB) = [];

ZeroLogicC = KeypointC.Location==0;
KeypointC.Location(ZeroLogicC) = [];
KeypointC.Level(ZeroLogicC) = [];
KeypointC.Scale(ZeroLogicC) = [];
KeypointC.Sign(ZeroLogicC) = [];

ZeroLogicD = KeypointD.Location==0;
KeypointD.Location(ZeroLogicD) = [];
KeypointD.Level(ZeroLogicD) = [];
KeypointD.Scale(ZeroLogicD) = [];
KeypointD.Sign(ZeroLogicD) = [];





end






% 
% 
% [maxval, maxind] = max(DoG,[],2);
% 
% [minval, minind] = min(DoG,[],2);
% 
% ZeroMax = (maxind <= 87) | (maxind == 2999);
% 
% ZeroMin = (maxind <= 87) | (maxind == 2999);
% 
% 
% LocationsMin = 1:PointCloud.LocationCount;
% LocationsMax = 1:PointCloud.LocationCount;
% LocationsMin(ZeroMin) = [];
% LocationsMax(ZeroMax) = [];
% 
% 
% MaxIndRemove = maxind(~ZeroMax);
% MaxValRemove = maxval(~ZeroMax);
% 
% MinIndRemove = minind(~ZeroMin);
% MinValRemove = minval(~ZeroMin);
% 
% 
% 
% for i = 1 : length(MaxIndRemove)
%     
%     CurrentNeighbors = Neighbors.Distance{LocationsMin(i)};
%     
% end
% 
% 
% 
% for i = 1 : length(MinIndRemove)
%     
%     
%     
% end
% 



















