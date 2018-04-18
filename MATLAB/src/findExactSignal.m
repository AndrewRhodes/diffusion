




function ExactSignal = findExactSignal(sigma, Theta, Phi, MaxDegreeL)

LengthTheta = length(Theta);
LengthPhi = length(Phi);

ExactSignal = zeros(length(Theta),1);




for iDegreeL = 0 : MaxDegreeL
    
    Coef = zeros( 2 * iDegreeL + 1, 1);
    Coefcsin = zeros( 2 * iDegreeL + 1, LengthTheta );
    
    for iOrderM = -iDegreeL : iDegreeL
        
        
        Coef(iOrderM + iDegreeL + 1 , 1) = (-1)^(iOrderM) * sqrt( ( (2*iDegreeL+1) / (4*pi) ) * ( factorial( iDegreeL - abs(iOrderM) ) / factorial( iDegreeL + abs(iOrderM) ) ) );
        
        
        if iOrderM < 0
            
            Coefcsin(iOrderM + iDegreeL + 1, :) = sqrt(2) * sin(abs(iOrderM) * Phi);
            
        elseif iOrderM > 0
            
            Coefcsin(iOrderM + iDegreeL + 1, :) = sqrt(2) * cos(iOrderM * Phi);
            
        else % iOrderM == 0
            
            Coefcsin(iOrderM + iDegreeL + 1, :) =  ones(LengthTheta,1);
            
        end

    end
    
        
    CoefExpDegreeL = exp( - (iDegreeL^2) / 9 );
    
    CoefExp = exp( - iDegreeL * (iDegreeL + 1) * (sigma^2)/2);
        
    Plm = legendre(iDegreeL, cos(Theta));
    
    Plm = [Plm(2:end,:); Plm];
    
    ExactSignal = ExactSignal + sum(CoefExpDegreeL * CoefExp * (Coef.*ones(2*iDegreeL+1, LengthTheta)) .* Plm .* Coefcsin ,1)';
    
end

ExactSignal = (20/(3*pi)) * ExactSignal;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for iDegreeL = 0 : MaxDegreeL - 1
%     
%    CoefExp = exp(-iDegreeL*(iDegreeL+1)*sigma^2/2); 
%     
%    ExactSignal = ExactSignal + sum(CoefExp * SignalOriginal{iDegreeL+1,1},1)';
%     
%     
% end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% YlmCoef = sqrt( ((2*DegreeL+1) / (4*pi)) * (( factorial( DegreeL - abs(0) ) / factorial( DegreeL + abs(0) ) )) );
% 
% Plm = legendre(DegreeL, cos(Theta));
% 
% Ylm = YlmCoef * Plm(1,:);
% 
% ExactSignal = Ylm * exp(-DegreeL*(DegreeL+1)*sigma^2/2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% for iDegreeL = 0 : MaxDegreeL - 1
%     
%     Coef = sqrt( (2*iDegreeL+1) / (4*pi) );
%     
%     CoefExp = exp( - iDegreeL * (iDegreeL + 1) * (sigma^2/2) );
%     
%     Plm = legendre(iDegreeL, cos(Theta));
%     
%     YlmCoef = sqrt( ((2*iDegreeL+1) / (4*pi)) * (( factorial( iDegreeL - abs(0) ) / factorial( iDegreeL + abs(0) ) )) );
% 
%     Ylm = YlmCoef * exp(1i*0*Phi) .* Plm(1,:)';
%     
%     
%     ExactSignal = ExactSignal + Coef * CoefExp * Ylm;
% end
% 
%     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% CoefSum = 0;
% 
% 
% for iDegreeL = 0 : MaxDegreeL - 1
%     
%     Coef = sqrt( (2*iDegreeL+1) / (4*pi) );
%     
%     CoefExp = exp( - iDegreeL * (iDegreeL + 1) * (sigma^2/2) );
% %         Coef*CoefExp
%     CoefSum = CoefSum + Coef*CoefExp;
%     
% end
% 
% ExactSignal = CoefSum * sum(cell2mat(OriginalSignal),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% for iDegreeL = 0 : MaxDegreeL - 1
%     
%     Coef = sqrt( (2*iDegreeL+1) / (4*pi) );
%     
%     CoefExp = exp( - iDegreeL * (iDegreeL + 1) * (sigma^2/2) );
%         
%     ExactSignal = ExactSignal + Coef*CoefExp*OriginalSignal{iDegreeL+1,1}';
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end