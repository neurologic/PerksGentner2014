function hyp = GOF(y,fitline)
%goodness of fit
hyp = 0;

(y - fitline)    % Errors
(y - fitline).^2   % Squared Error
mean((y - fitline).^2)   % Mean Squared Error
RMSE = sqrt(mean((y - fitline).^2));  % Root Mean Squared Error


sampsize = max(size(y));
for irand = 1:20  %get distribution of random estimates for RMSE
randint = randsample(sampsize,sampsize);
randy = y(randint);
RMSErand(irand) = sqrt(mean((randy - fitline).^2));  
end

%is the RMSE from the data different from the estimate off of rand data
% --> is it outside of 2SD from mean?

q = quantile(RMSErand,[0.15,0.85]);

if RMSE < q(1) % if the RMSE is less than the 15% ci on the RMSErand 
    
hyp = 1;
end

if RMSE > q(2)
%(lower RMSE is better fit... so the fit COULD be worse than the rand fit)
    
hyp = -1;
end