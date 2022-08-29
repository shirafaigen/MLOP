%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the Manifold Locally Optimal Projection (MLOP) d for denoising and reconstruction of low-dimensional
% manifold in high-dimensional space. We suggest a multidimensional extension of the Locally
% Optimal Projection algorithm which was introduced by Lipman et al. in 2007 for
% surface reconstruction in 3D. The method bypasses the curse of dimensionality and avoids
% the need for carrying out dimensional reduction. It is based on a non-convex optimization
% problem, which leverages a generalization of the outlier robust L1-median to higher
% dimensions while generating noise-free quasi-uniformly distributed points reconstructing
% the unknown low-dimensional manifold. 

% The MLOP framework generate a manifold object, that inherits from the
% Data_Interface (e.g. Data_Cylinder_Non_Uniform_Sample). 
% In case a new manifold is used, the new class should inherit from "Data_Interface", and implement the methods:
%  1) getInitialH
%  2) createData2
%  3) createData_with_resolution
%  4) drawData

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MLOP algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The MLOP follows the following steps (for additional information please see [1]):
%   Input:  $P=\{p_j\}_{j=1}^J \subset \mathbb{R}^n$, $\epsilon>0$
%   Output: $Q=\{q_i\}_{i=1}^I \subset \mathbb{R}^n$

%   Initialize $Q^{(0)}$ as a subsample of $P$
%   Estimate $h_1$ and $h_2$
%   Repeat
%       For each $q_{i'}^{(k)} \in Q^{(k)}$
%           Calculate $\nabla G(q_{i'}^{(k)})$ by assessing $\alpha_j^{i'}$, $\beta_i^{i'}$}
%           q_{i'}^{(k+1)} = q_{i'}^{(k)} - \gamma_k \nabla G(q_{i'}^{(k)})$}
%       EndFor
%   Until{$\|\nabla G(q_{i'}^{(k)})\| <\epsilon$ or max number of iterstions is reached}

% The end product of the function is the set q - that reconstructs the manifold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   References:
%      [1] Lipman, Y., Cohen-Or, D., Levin, D., Tal-Ezer, H.
%          Parameterization-free projection for geometry reconstruction.
%          In: ACM Transactions on Graphics (TOG), vol. 26, p. 22. ACM (2007)

%      [2] Faigenbaum-Golovin, S., & Levin, D. (2020). 
%          Manifold Reconstruction and Denoising from Scattered Data in High Dimension via a Generalization of L1-Median. arXiv preprint arXiv:2012.12546.?

% See also MLOP enhancments
%      [3] Faigenbaum-Golovin, S., & Levin, D. (2020). 
%           Approximation of Functions over Manifolds in High Dimension from Noisy Scattered Data . arXiv preprint arXiv:2012.13804.?

%      [4] Faigenbaum-Golovin, S., & Levin, D. (2020). 
%           Manifold Repairing, Reconstruction and Denoising from Scattered Data in High-Dimension. arXiv preprint arXiv:2012.13804.?

%   (C) Faigenbaum-Golovin, S., & Levin, D. Tel-Aviv University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear df_k;
clear w
clear w1
clear w2
clear w3

close all;

GENERATE_DATA = false; %Generate the manifold data, the result will be stored in a vector P, or use the previosly generated data

global d
d = 0.1; % the regularisation parameter fot eh Hd norm calculation
debug_index = 1; % index for saving debugging information
extrinsic_dim = 60; % the extrinsic dimension of the data
intrinsic_dim = -1; % the intrinsic dimension of the data - not used
approx_point_num_p = 1000; % approximate number of how many points will be generated
add_noise = false; %add noise along with the data generation, or add it later 
num_points_in_q = 200; % number of points in the reconstruction set Q

UP_SAMPLE_RECONSTRUCTION = false;
GD_num_iter = 1000; % number of iterations to perform the MLOP gradient descent (a max value)
DIM_RED = true; % should we generate the linear sketching Dim. reduction matrix, or calculate norm in the current dimension by using the eye matrix
new_dim = 20; % the dimension, where we project the data to.
ESTIMATE_H = true; % extimate the fill distance?

% Generate the manifold data, the result will be stored in a vector P
% it is possible to generate the data once, and disable this block for the
% following executions
if GENERATE_DATA

    %%%%%%%%%%%%%%  Create the Cylinder object and Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     data_obj = Data_Cylinder_Non_Uniform_Sample();

%         data_obj = Data_Cylinder_Non_Uniform_Sample();
%         data_obj = Data_AncientGeneratedLetters();
%         data_obj = Data_Cone_Non_Uniform_Sample(); h1 = 0.35; h2 = h1;
%           data_obj = Data_Orthogonal_Matrix_Uniform_Sample(); 
%         data_obj = Data_Elipse(); add_noise_to_sampled_data = false;
%         data_obj = Data_RN_Cylinder_Non_Uniform_Sample();
%         data_obj = Data_Plane_Uniform_Sample();
%         data_obj = Data_Functions_Sample();



    p = data_obj.createData(extrinsic_dim, intrinsic_dim, approx_point_num_p, add_noise);
	num_points_in_p = size(p,1);
    
	%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	pp_noiseless = p;
	ones_p = ones(size(p, 1), extrinsic_dim);
	noise = 0.1; %0.2;
   
	p = pp_noiseless + noise*(rand(num_points_in_p, extrinsic_dim) - 0.5*ones_p); 
     
	%%%%%%%%%%%%%% Randomly sampe the set Q for the reconstruction %%%%%%%%%%%%%%%%%

	q_1 = p(randsample(num_points_in_p, num_points_in_q, false), :);   
    q_1 = unique(q_1,'rows');
    num_points_in_q = size(q_1,1);
    
     % Upsample the point-set P, with replacment, and add a bit of noise each time
    if (UP_SAMPLE_RECONSTRUCTION)

        p1 = ones(num_points_in_q, extrinsic_dim);
        qq_1 = q_1;
        num_of_duplicates =1;
        for i=1:num_of_duplicates%10
            noise = 0.2;
            qq_t = q_1 + noise*(rand(num_points_in_q, extrinsic_dim) - 0.5*p1); % 315x4
            qq_1 = [qq_1;qq_t];
        end
        q_1 = qq_1;
    end

	%%%%%%%%%%%%%% Generate the linear sketching Dim. reduction matrix %%%%%%%%%%%%%%%%%

    if (DIM_RED)
        % perform linear global sketching, by projecting the data on a
        % normal distribution matrix, and performing QR decomposition of it.
        mu = 0;
        sigma = 1;
        
        pp_t = p';
        G = normrnd(mu,sigma,[size(pp_t, 2), new_dim]);
        B = pp_t*G;
        [Q,R] = qr(B, 0);
        pp_n = p*Q;
        DimRedM = Q;
    else % or use the high dimensional data as is
        DimRedM = eye(size(p, 2));
    end

    % Visualize the data
    data_obj.drawData(p, q_1, 10, intrinsic_dim)
     
    %%%%%%%%%%%%%% Estimate the fill-distances of P and Q (h1, h2 respectivly) %%%%%%%%%%%%%%%%% 
    Sp = round(size(p, 1)/size(q_1, 1)); % extimate the number of Ps that are served with a q_i

    if (ESTIMATE_H)
        
        minNumOfPoints = 8*Sp;
        %%%%%%%%%%%%%% Estimate fill-distance of P %%%%%%%%%%%%%%
        % get the initial guess for h, where we should start looking 
        initial_h_1  = data_obj.getInitialH(p, DimRedM, minNumOfPoints); 
        h1 = calculateFillDistancePerPoint(p, DimRedM, minNumOfPoints, true, initial_h_1, []);

        h1 = median(h1);
        if (sum(isnan(h1)) > 0)
            return;
        end

        %%%%%%%%%%%%%% Estimate fill-distance of Q %%%%%%%%%%%%%%
        % randomly and uniformly sample a new set Q, that will represent
        % the optimal Q*, and calculate its fill-disatnce
        index = randsample(num_points_in_p, num_points_in_q);
        index = sort(index);
        qq_uniform = p(index,:);
        initial_h_2  = data_obj.getInitialH(qq_uniform, DimRedM, minNumOfPoints/2); 
        
        h2 = calculateFillDistancePerPoint(qq_uniform, DimRedM, minNumOfPoints/2, false, initial_h_2, []);
        h2 = median(h2);
    end
end

q = q_1;
       
% Visualize the data
data_obj.drawData(p, q, 10, intrinsic_dim);
debug_index_2 = 1;

%%%%%%%%%%%%%% Start MLOP reconstruction iterations %%%%%%%%%%%%%%%%% 
for j=1:GD_num_iter
    
    % Once in a while write the iteration index
    if (mod(j, 50) == 0)
        j
	end
      
    % In case we are not in the first iteration, update the q_k, and df_k
    % parameters for the GD step size calculation
	if (j>1)
        q_k_1 = q_k;
        df_k_1 = df_k;
    end
   
	q_k = q;
     
    %%%%%%%%%%%%%%%%% For each q_i, calculate its GD step %%%%%%%%%%%%%%%%%
    for i=1:size(q,1)
        
        %%%%%%%%%%%%%%%%% Calculate the distances of the points in Q, and P from the current q_i %%%%%%%%%%%%%%%%%
        odp = repmat(q(i,:), size(p, 1), 1);     
        odq = repmat(q(i,:), size(q, 1), 1);

        eeq = calculateNorm(q(i, :), q, DimRedM);
        eep = calculateNorm(q(i, :), p, DimRedM);


        %%%%%%%%%%%%%%%%% Calculate the derivatives %%%%%%%%%%%%%%%%%
         Hd_p = Hd(eep,d);
        hdvec = (theta(eep, h1)./Hd_p).*(1 - 2./h1^2 * Hd_p.^2); 

        [eta, d_eta] = eta_function(eeq);

        betvec = 2*(theta(eeq, h2)./eeq).*(abs(d_eta) + 1/h2^2 *eta.*eeq ); % redivative by varios eta, derive with respect to q,
        betvec(i,1)=0;
        
        v = hdvec'*(odp-p)/sum(hdvec); % the E_1 term derivative
        u = betvec'*(odq-q);           % the E_2 term derivative

        %%%%%%%%%%%%%%%% In case we are in the first iteration, calculate the balancing term (lambda) %%%%%%%%%%%%%%%%
        if (j==1)
            norm_v = sqrt(v*v');
            norm_u = sqrt(u*u');
            
            w1(i) = 1/norm_v; 
            w2(i) = 2*1/norm_u;
        end

        %%%%%%%%%%%%%%%% Calculate the gradient %%%%%%%%%%%%%%%%
        df_k(i,:) =  w1(i)*v  - w2(i)*u ; 
        
%         df_k(i,:) = 1/(sum(hdvec))*(v - w(i)*u);% 1/(sum(hdvec))* SHIRA


        %%%%%%%%%%%%%%%% Update the gamma Gradient Descent coeficient    %%%%%%%%%%%%%%%%%%%%%%
        if (j<2)
            gama(i) = 0.1; % an initial gamma
        else
            d_q = q(i,:) - q_k_1(i,:);
            d_deriv = df_k(i,:) - df_k_1(i,:);
        
            gama(i) = (d_q*d_deriv')/(d_deriv*d_deriv');
            if ( isnan(gama(i)))
                 gama(i) = 0;
                 ['const gamma ' num2str(i) ' ' num2str(j)]
            end
        end
        
        %%%%%%%%%%%%%%%% Calc GD iteration as: q^(k+1) = q^(k) + gamma*df/dq(E1+E2) %%%%%%%%%%%
        next_q_i= q(i,:) - gama(i)*df_k(i,:);
        
        
        %%%%%%%%%%%%%%%% In case the q_i is nan, indicate this to the uset %%%%%%%%%%%%%%%%
        nanq = isnan(next_q_i);
        if (sum(nanq(:)) > 0)
           'Error'
           
           % Visualize the data
            data_obj.drawData(p, q, 9, intrinsic_dim);
        end
       
        q(i,:) = next_q_i;

        % Collect some statistic about the number of points in the support, for debugging
        med_data_min_points_q(debug_index_2, 1:5) = [w1(i)*median(v), w2(i)*median(u), sum(abs(hdvec)>0.01), sum(abs(betvec)> 0.01),  1];
        debug_index_2 = debug_index_2 + 1;
    end
    
    debug_index_2 = 1;
    
    %%%%%%%%%%%%%%%% In the first iteration update the coef %%%%%%%%%%%%%%%%
    if (j==1)
        w1 = median(w1)*ones(size(w1));
        w2 = median(w2)*ones(size(w2));
    end
    
    % Collect some statistic about the number of points in the support, for debugging
    med_data_min_points(debug_index, 1:5) = median(med_data_min_points_q(:, 1:5));
    debug_index = debug_index+1;
    
    if (mod(j, 10) == 0)
        
        % Visualize the data
        data_obj.drawData(p, q, 11, intrinsic_dim);


        % [hd, relative_err_mean, relative_err_max] =
        % data_obj.evaluateError(p, q, DimRedM); % calculate the error
        'l';
    end

end

'MLOP finished'