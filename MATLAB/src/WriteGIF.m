function WriteGIF(FileName, i, T)
%% Generate GIF from frames
%
% WRITEGIF produces a GIF from the current frame. A new GIF file is created
% upon calling WRITEGIF initially, and each subsequent step appends the new
% frame to the .gif file.
%
%
% Input: 
%       FileName     =  String containing the desired filename
%       i            =  Frame number (must start at 1)
%       T            =  Delay time between frames
%
%
% Output: 
%       None
%
%
% Do not distribute without permission from AFRL/RVSV GNC Group
%
%

%% Convert current frame to image
% M = getframe(gcf,[0,0,1850,1050]);
ax = gca;
ax.Units = 'pixel';
pos = ax.Position;
ti = ax.TightInset;
rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
M = getframe(ax, rect);
im=frame2im(M);
% im = imresize(im, 0.6);

[imind,cm] = rgb2ind(im,256);


%% Write GIF
if i == 1;
    imwrite(imind,cm, FileName, 'gif', 'Loopcount', inf, 'DelayTime', T);
else
    imwrite(imind, cm, FileName, 'gif', 'WriteMode', 'append', 'DelayTime', T);
end


end