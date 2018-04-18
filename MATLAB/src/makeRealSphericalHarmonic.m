




function SphericalHarmonic = makeRealSphericalHarmonic( MaxDegreeL, Theta, Phi )


LengthTheta = length(Theta);
LengthPhi = length(Phi);
if LengthTheta ~= LengthPhi
    error('Theta and Phi must have same dimensions')
end

Ylm = cell( MaxDegreeL, 1 );
% Ylm = zeros(1, LengthTheta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iDegreeL = 0 : MaxDegreeL
    
    Coef = sqrt( ( (2*iDegreeL+1) / (4*pi) ) );
    
    Plm = legendre(iDegreeL, cos(Theta));
    
    Plm = Plm(1,:);
    
    Ylm{iDegreeL + 1, 1} = Coef .* Plm;
    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for iDegreeL = 0 : MaxDegreeL - 1
%         
%     
%     Coef = sqrt( ( (2*iDegreeL+1) / (4*pi) ) );
%        
%     Plm = legendre(iDegreeL, cos(Theta));
%     
%     Plm = Plm(1,:);
%     
%     Ylm{iDegreeL + 1, 1} = Coef .* Plm;
%     
% end
% 
% 
% SphericalHarmonic = Ylm;








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for iDegreeL = 0 : MaxDegreeL - 1
%     
%     Coef = zeros( 2 * iDegreeL + 1, 1);
%     Coefcsin = zeros( 2 * iDegreeL + 1, LengthTheta );
%     
%     for iOrderM = -iDegreeL : iDegreeL
%         
%         
%         Coef(iOrderM + iDegreeL + 1 , 1) = (-1)^(iOrderM) * sqrt( ( (2*iDegreeL+1) / (4*pi) ) * ( factorial( iDegreeL - abs(iOrderM) ) / factorial( iDegreeL + abs(iOrderM) ) ) );
%         
%         
%         if iOrderM < 0
%             
%             Coefcsin(iOrderM + iDegreeL + 1, :) = sqrt(2) * sin(abs(iOrderM) * Phi);
%             
%         elseif iOrderM > 0
%             
%             Coefcsin(iOrderM + iDegreeL + 1, :) = sqrt(2) * cos(iOrderM * Phi);
%             
%         else % iOrderM == 0
%             
%             Coefcsin(iOrderM + iDegreeL + 1, :) =  ones(LengthTheta,1);
%             
%         end
% 
%     end
%     
%     
%     Plm = legendre(iDegreeL, cos(Theta));
%     
%     Plm = [Plm(2:end,:); Plm];
%     
%     Ylm{iDegreeL + 1, 1} = (Coef.*ones(2*iDegreeL+1, LengthTheta)) .* Plm .* Coefcsin ;
%     
%     
% end
% 
% 
% SphericalHarmonic = Ylm;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% for iDegreeL = 0 : MaxDegreeL - 1
%     
%     Coef = zeros( 2 * iDegreeL + 1, 1);
%     CoefExp = zeros( 2 * iDegreeL + 1, LengthTheta );
%     
%     for iOrderM = -iDegreeL : iDegreeL
%         
%         Coef(iOrderM + iDegreeL + 1 , 1) = (-1)^((iOrderM + abs(iOrderM))/2) * sqrt( ( (2*iDegreeL+1) / (4*pi) ) * ( factorial( iDegreeL - abs(iOrderM) ) / factorial( iDegreeL + abs(iOrderM) ) ) );
%         CoefExp(iOrderM + iDegreeL + 1, :) = exp( 1i * iOrderM * Phi'); % Real & Imag      
%         
%     end
%     
%     Plm = legendre(iDegreeL, cos(Theta));
%     
%     Plm = [Plm(2:end,:); Plm];
%     
%     Ylm{iDegreeL + 1, 1} = (Coef.*ones(2*iDegreeL+1, LengthTheta)) .* Plm .* CoefExp ;
%     
% %     pause
% end
% 
% SphericalHarmonic = cellfun(@real, Ylm, 'UniformOutput', 0);
% % SphericalHarmonic = Ylm;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for iDegreeL = 0 : MaxDegreeL - 1
%     
%     Coef = zeros(iDegreeL+1, 1);
%     CoefExp = zeros(iDegreeL+1,LengthTheta);
% %     CoefCos = zeros(iDegreeL+1,LengthTheta);
%     
%     for iOrderM = 0 : iDegreeL
%         
%         Coef(iOrderM + 1,1) = sqrt( ( (2*iDegreeL+1) / (4*pi) ) * ( factorial( iDegreeL - iOrderM ) / factorial( iDegreeL + iOrderM ) ) );
%         
% %         CoefCos(iOrderM + 1,:) = cos(iOrderM * Phi); % Real
%         CoefExp(iOrderM + 1,:) = exp( 1i * iOrderM * Phi); % Real & Imag
%         
%     end
%     
%     Plm = legendre(iDegreeL, cos(Theta));
%     Ylm{iDegreeL+1, 1} = (Coef.*ones(iDegreeL+1, LengthTheta)) .* Plm .* CoefExp ;
% %     Ylm{iDegreeL+1, 1} = sqrt(2) * (Coef.*ones(iDegreeL+1, LengthTheta)) .* Plm .* CoefCos ;
%     
% end
% 
% SphericalHarmonic = cellfun(@real, Ylm, 'UniformOutput', 0);
% % SphericalHarmonic = Ylm;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






end