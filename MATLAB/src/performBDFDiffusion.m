% Andrew Rhodes
%
% input[SignalIn]: the signal to be diffused. size[m x 1]
% input[NumSteps]: The number of steps of diffusion. scalar
% input[ItL]: (I - a * t * L) size[m x m]
% output[SignalOut]: the diffused signal over all time steps. 
%                    size[m x NumSteps]



function SignalOut = performBDFDiffusion(SignalIn, NumSteps, ItL)

NumVertex = length(SignalIn);

OrderBDF = length(ItL);


SignalOut = zeros(NumVertex, NumSteps);
SignalOut(:,1) = SignalIn;


if OrderBDF == 1
    SignalOut = runBDF1(SignalOut, NumSteps, ItL{1,1});
elseif OrderBDF == 2
    SignalOut = runBDF2(SignalOut, NumSteps, ItL{1,1}, ItL{2,1});
elseif OrderBDF == 3
    SignalOut = runBDF3(SignalOut, NumSteps, ItL{1,1}, ItL{2,1}, ItL{3,1});
elseif OrderBDF == 4
    SignalOut = runBDF4(SignalOut, NumSteps, ItL{1,1}, ItL{2,1}, ItL{3,1}, ItL{4,1});
else
    error('Expected OrderBDF = 1, 2, 3, 4, but received %i', OrderBDF)
end





end


function SignalOut = runBDF1(SignalOut, NumSteps, ItL)


%WaitBar = waitbar(0, sprintf('BDF1 Diffusion %i of %i', 0, NumSteps - 1));

for i = 1 : NumSteps - 1
    
    
    [SignalOut(:,i+1), flag] = bicg(ItL, SignalOut(:,i), 1e-10, 100);
    
    
    if flag
        disp(flag)
    end
    
%	  waitbar(i/NumSteps, WaitBar, sprintf('BDF1 Diffusion %i of %i', i, NumSteps-1));
end

%waitbar(i/NumSteps, WaitBar, sprintf('BDF1 Diffusion Complete'));
%close(WaitBar)

end




function SignalOut = runBDF2(SignalOut, NumSteps, ItL, I23tL)


%WaitBar = waitbar(0, sprintf('BDF2 Diffusion %i of %i', 0, NumSteps - 1));

for i = 1 : NumSteps - 1
    
    if i == 1
        [SignalOut(:,i+1), flag] = bicg(ItL, SignalOut(:,i), 1e-10, 100);
    else
        [SignalOut(:,i+1), flag] = bicg(I23tL, (4/3)*SignalOut(:,i) - (1/3)*SignalOut(:,i-1), 1e-10, 100);
    end
    
    if flag
        disp(flag)
    end
    
%     waitbar(i/NumSteps, WaitBar, sprintf('BDF2 Diffusion %i of %i', i, NumSteps-1));
end

%waitbar(i/NumSteps, WaitBar, sprintf('BDF2 Diffusion Complete'));
%close(WaitBar)


end




function SignalOut = runBDF3(SignalOut, NumSteps, ItL, I23tL, I611tL)

% WaitBar = waitbar(0, sprintf('BDF3 Diffusion %i of %i', 0, NumSteps - 1));

for i = 1 : NumSteps - 1
    
    if i == 1
        [SignalOut(:,i+1), flag] = bicg(ItL, SignalOut(:,i), 1e-10, 100);
    elseif i == 2
        [SignalOut(:,i+1), flag] = bicg(I23tL, (4/3)*SignalOut(:,i) - (1/3)*SignalOut(:,i-1), 1e-10, 100);
    else
        [SignalOut(:,i+1), flag] = bicg(I611tL, (18/11)*SignalOut(:,i) - (9/11)*SignalOut(:,i-1) + (2/11)*SignalOut(:,i-2), 1e-10, 100);
    end
    
    if flag
        disp(flag)
    end
    
%     waitbar(i/NumSteps, WaitBar, sprintf('BDF3 Diffusion %i of %i', i, NumSteps-1));
end

% waitbar(i/NumSteps, WaitBar, sprintf('BDF3 Diffusion Complete'));
% close(WaitBar)

end





function SignalOut = runBDF4(SignalOut, NumSteps, ItL, I23tL, I611tL, I1225tL)

% WaitBar = waitbar(0, sprintf('BDF4 Diffusion %i of %i', 0, NumSteps - 1));

for i = 1 : NumSteps - 1
    
    if i == 1
        [SignalOut(:,i+1), flag] = bicg(ItL, SignalOut(:,i), 1e-10, 100);
    elseif i == 2
        [SignalOut(:,i+1), flag] = bicg(I23tL, (4/3)*SignalOut(:,i) - (1/3)*SignalOut(:,i-1), 1e-10, 100);
    elseif i == 3
        [SignalOut(:,i+1), flag] = bicg(I611tL, (18/11)*SignalOut(:,i) - (9/11)*SignalOut(:,i-1) + (2/11)*SignalOut(:,i-2), 1e-10, 100);
    else
        [SignalOut(:,i+1), flag] = bicg(I1225tL, (48/25)*SignalOut(:,i) - (36/25)*SignalOut(:,i-1) + (16/25)*SignalOut(:,i-2) - (3/25)*SignalOut(:,i-3), 1e-10, 100);
    end
    
    if flag
        disp(flag)
    end
    
%     waitbar(i/NumSteps, WaitBar, sprintf('BDF4 Diffusion %i of %i', i, NumSteps-1));
end

% waitbar(i/NumSteps, WaitBar, sprintf('BDF4 Diffusion Complete'));
% close(WaitBar)


end




















