

function D_hat = findOptimalDoG(J, H, c, xhat)

D_hat = c + J' * xhat + xhat' * (H \ xhat);


end



