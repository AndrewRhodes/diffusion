

function [J, H, c] = buildJHc(y)

J = y(1:4,1);

H = [y(5), y(6),  y(7),  y(8);
    y(6), y(9),  y(10), y(11);
    y(7), y(10), y(12), y(13);
    y(8), y(11), y(13), y(14)];

c = y(15);

end