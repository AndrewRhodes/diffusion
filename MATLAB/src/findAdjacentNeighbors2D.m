



function Neighbors = findAdjacentNeighbors2D(ImageSize)

NumPixels = prod(ImageSize);



[Iindex, Jindex] = ind2sub([ImageSize(1),ImageSize(2)], (1:NumPixels)');

Neighbors = cell(NumPixels,1);


for i = 1 : NumPixels
    
%     CurrentPoint = [Iindex(i), Jindex(i)];
    
    PossibleNeighbors = [Iindex(i)+1, Jindex(i);
                        Iindex(i)-1, Jindex(i);
                        Iindex(i), Jindex(i)+1;
                        Iindex(i), Jindex(i)-1;
                        Iindex(i)+1, Jindex(i)+1;
                        Iindex(i)+1, Jindex(i)-1;
                        Iindex(i)-1, Jindex(i)+1;
                        Iindex(i)-1, Jindex(i)-1];
                    
                    
    % Remove neighbors from outside image bounds
    PossibleNeighbors(PossibleNeighbors(:,1) > ImageSize(1), :) = [];
    PossibleNeighbors(PossibleNeighbors(:,1) < 1, :) = [];
    PossibleNeighbors(PossibleNeighbors(:,2) > ImageSize(2), :) = [];
    PossibleNeighbors(PossibleNeighbors(:,2) < 1, :) = [];
    
    
    Neighbors{i,1} = PossibleNeighbors;
    
end




end