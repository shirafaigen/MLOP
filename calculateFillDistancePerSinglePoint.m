%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function estimate the fill-distance for the "currentPoint" given with
% respect to the P set.
% The function go overt the h_s set, and finds the minimal h, that will
% give at least #minNumOfPoints points in the "support" of the Gaussian
% with respect to the E_1, E_2 term (as indicated by the byFirstTerm
% parameter).

% The function returns the optimal h for the current point with respect to
% P, or NAN, if no h was found that gives enough points in the support

% Input parameters:
% p - the set of noisy data ref data, 
% currentPoint - current point for fill-distance calculation
% DimRedM - dimensional reduction matrix for the linear sketching
% minNumOfPoints - lower bound on the min number of points in the support of the Gausisan function
% byFirstTerm - true/false, optimize h by E_1 term or E_2
% h_s - a list of possiable fill-distances to test
% i - the index of the current point
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h, num_points_in_support] = calculateFillDistancePerSinglePoint(p, currentPoint, DimRedM, minNumOfPoints, byFirstTerm, h_s, i)
global ALPHA_BETA_FROM_PAPER;
global d;
h = nan;
num_points_in_support = nan;
for j=1:size(h_s, 2)
         
    if (~isnan(h))
        continue;
    end

    h1 = h_s(j);
    h2 = h_s(j);

    % calculate the norms, and the gradient of the first/second term.
    % check that in the support of this term, there are enough points
    eep = calculateNorm(currentPoint, p, DimRedM); 
    eeq = calculateNorm(currentPoint, p, DimRedM);

        
    if (ALPHA_BETA_FROM_PAPER)
        alphavec = theta(eep,h1)./eep; % from the paper
        alphavec(alphavec == Inf) = 0; %note, that if q is randsample from p. then p(i)~=q(i) and eqq(666) = 0 and not in the i's place


        betvec=(theta(eeq,h2)./eeq).*(1./(eeq/h2).^4);  % from the paper       
        betvec(i,1)=0;

    else

        Hd_p = Hd(eep,d);
        alphavec = (theta(eep, h1)./Hd_p).*(1 - 2./h1^2 * Hd_p.^2); % the real derive with respect to q,c and Hd

        [eta, d_eta] = eta_function(eeq);

        betvec = (theta(eeq, h2)./eeq) .* (abs(d_eta) + 2/h2^2 *eta.*eeq ); % redivative by varios eta, derive with respect to q,
        betvec(i,1)=0;

    end
        
    if (byFirstTerm)
        
        if (sum(abs(alphavec)>0.01) > minNumOfPoints ) % there shold be at least "5" points for the approximation
            h = h1;
            num_points_in_support = sum(abs(alphavec)>0.01);
        end
    else

        if (sum(abs(betvec)> 0.01) > minNumOfPoints ) % there shold be at least "5" points for the approximation
            h = h1;
            num_points_in_support = sum(abs(betvec)>0.01);
        end
    end
end
'l';