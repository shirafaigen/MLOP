%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function estimate the fill-distance of the given set P.
% The function uses the "initial_h_guess" parameter, to gemerate a neighborhood of fill-distances (h_s).
% later for each points in P, it go overt this h_s set, and finds the minimal h, that will
% give at least #minNumOfPoints points in the "support" of the Gaussian
% with respect to the E_1, E_2 term (as indicated by the byFirstTerm parameter).

% The function returns the optimal h for the current point with respect to
% P, or NAN, if no h was found that gives enough points in the support

% Input parameters:
% p - the set of noisy data ref data, 
% DimRedM - dimensional reduction matrix for the linear sketching
% minNumOfPoints - lower bound on the min number of points in the support of the Gausisan function
% byFirstTerm - true/false, optimize h by E_1 term or E_2
% initial_h_guess - the initial guess for the fill-distance, which aids in
% calculating different h_s in its surrounding
% currentPoint - a reference points for which we calculate the h. id the
% points is empty, the fill-distance is estimated for the entier set P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = calculateFillDistancePerPoint(p, DimRedM, minNumOfPoints, byFirstTerm, initial_h_guess, currentPoint) 
N = size(p, 1);

% Use the initial guess of the h function, to create a list of h_s values
% (in the length of 40)  in the range of [h_0*(1-epsilon), h_0*(1+epsilon)], with the stepsize of "step"
per_p = 0.8; 
min_h = max(0.01, initial_h_guess*(1-per_p));
max_h = initial_h_guess*(1+per_p);
step = (max_h - min_h)/40;

% make this list global, to use in the inner function 
h_s = min_h:step:max_h;

h_candidates = [];
['Total amount of hs to eval:' num2str(size(h_s, 2))]
% if there is a reference point, calculate the fill-distance with respect
% to it, otherwise for the entier set P
if (~isempty(currentPoint))
    
    h = calculateFillDistancePerSinglePoint(p, currentPoint, DimRedM, minNumOfPoints, byFirstTerm, h_s, 1);
    
    % indicate if a suitable fill-distance was not found
    if (isnan(h))
        'Error there is a NAN h'
        
        h = initial_h_guess;
    end
    
else

    % find the representative fill-distance for 20% of the data, for fact calculations
    num_to_sample = round(N*0.20);
    indexes = (round(rand(round(num_to_sample), 1)*(N-1))+1);
    
    for kk=1:num_to_sample
        i = indexes(kk);
       
        if (mod(kk, 10) == 0)
             kk
        end

        % for each point estimate the fill-distance, with respect to the
        % ones in the list h_s
        currentPoint = p(i,:);
        [h(kk), num_points_in_support(kk)] = calculateFillDistancePerSinglePoint(p, currentPoint, DimRedM, minNumOfPoints, byFirstTerm, h_s, i);
        if (isnan(h(kk)))
            'Error there is a NAN h'
            kk
            median([h(~isnan(h))])
            break;
        end
    end

end

% indicate what was the problem with finding the right h
if (sum(isnan(h)) > 0)
    'Error there is a NAN h'
    h
    median([h(~isnan(h))])
    h=nan;
end

if (sum(h == h_s(1)) > 0.2*size(h, 2))
    'Error initial h is too large'
    median([h(~isnan(h))])
    h = nan;
end

'l';