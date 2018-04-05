




function SphericalHarmonic = makeRealSphericalHarmonic( MaxDegreeL, Theta, Phi )


LengthTheta = length(Theta);
LengthPhi = length(Phi);
if LengthTheta ~= LengthPhi
    error('Theta and Phi must have same dimensions')
end

Ylm = cell( MaxDegreeL, 1 );
% Ylm = zeros(1, LengthTheta);

for iDegreeL = 0 : MaxDegreeL-1
    
    Coef = zeros(iDegreeL+1, 1);
    CoefExp = zeros(iDegreeL+1,LengthTheta);
%     CoefCos = zeros(iDegreeL+1,LengthTheta);
    
    for iOrderM = 0 : iDegreeL
        
        Coef(iOrderM + 1,1) = sqrt( ( (2*iDegreeL+1) / (4*pi) ) * ( factorial( iDegreeL - iOrderM ) / factorial( iDegreeL + iOrderM ) ) );
        
%         CoefCos(iOrderM + 1,:) = cos(iOrderM * Phi); % Real
        CoefExp(iOrderM + 1,:) = exp( 1i * iOrderM * Phi); % Real & Imag
        
    end
    
    Plm = legendre(iDegreeL, cos(Theta));
    Ylm{iDegreeL+1, 1} = (Coef.*ones(iDegreeL+1, LengthTheta)) .* Plm .* CoefExp ;
%     Ylm{iDegreeL+1, 1} = sqrt(2) * (Coef.*ones(iDegreeL+1, LengthTheta)) .* Plm .* CoefCos ;
    
end

SphericalHarmonic = cellfun(@real, Ylm, 'UniformOutput', 0);
% SphericalHarmonic = Ylm;

end