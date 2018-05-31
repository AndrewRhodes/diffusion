

function Phi = buildPhi(vx, vy, vz, sigma)

Phi = [vx^2, 2*vx*vy, 2*vx*vz, 2*vx*sigma, ...
    vy^2, 2*vy*vz, 2*vy*sigma, ...
    vz^2, 2*vz*sigma, ...
    sigma^2];


end