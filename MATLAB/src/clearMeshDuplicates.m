function [VerticesOut, FacesOut] = clearMeshDuplicates(VerticesIn, FacesIn)
% ----------------------------------------------------------------------------
% FUNCTION clearMeshDuplicates - Cleans up duplicate Mesh and faces.
%
% INPUTS:
%    	 INPUT1 - Vertices array 		(m1x3)
%    	 INPUT2 - Faces array			(m2x3, m1=m2*3)
% OUTPUTS:
%    	OUTPUT1 - New vertices array   	(m3x3)
%    	OUTPUT2 - New faces array		(m4x3, m3=m4*3)
% ----------------------------------------------------------------------------
% Created:     11/18/2015
% Author:      Eric Kim, Andrew Rhodes
% Department:  WVU Applied Space Exploration Laboratory (ASEL) 
%              West Virginia Robotic Technology Center
% Contact:     eric.kim@mail.wvu.edu, arhodes5@mix.wvu.edu
% Copyright:   2015-2016
% ----------------------------------------------------------------------------

% load sample_firstloop.mat
fprintf('Clearing out duplicate vertices and faces... \n'); tic

% Transform FacesIn into 1 column
FacesInT = FacesIn';
FacesInRow = FacesInT(:);
FacesLength = length(FacesInRow);

% Initialize Outputs
FaceArray = zeros(FacesLength, 1);
VerticesOut = zeros(FacesLength, 3);

% Counters
i1 = 1; faceN = 0;
while (i1 <= FacesLength)
    V_i1 = VerticesIn(FacesInRow(i1), :);
    % Find any rows after the index that equals the index vertices
    if (FaceArray(i1) == 0)
        VerticesOut(i1,:) = V_i1;
        faceN = faceN + 1;
        FaceArray(i1) = faceN;
        for i2 = i1 + 1 : FacesLength 
            if (VerticesIn(FacesInRow(i2), 1) == V_i1(1))
                if (VerticesIn(FacesInRow(i2), 2) == V_i1(2))
                    if (VerticesIn(FacesInRow(i2), 3) == V_i1(3))
                        FaceArray(i2) = faceN;
                    end
                end
            end
        end
    end
    i1 = i1 + 1;
end
fprintf('Completed! '); toc

%Transform back to original mx3 array
FacesOut = vec2mat(FaceArray, 3);
VerticesOut( ~any(VerticesOut, 2), :) = [];
end