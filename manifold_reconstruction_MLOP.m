GENERATE_DATA = true;
sp = [];
param_q = false;

clear df_k;
clear w
clear w1
clear w2
clear w3
clear w4
close all;
global d
d = 0.2;
kk1=1;


if GENERATE_DATA
%     p = create_manifold_repairing_data(); % returns normalised data to [0,1]
     data_obj = Data_Cylinder_Non_Uniform_Sample();
     p = data_obj.createData(60, 60, 1000, false);
%      p=p(:,1:3);
     
%      peak_point = [0.06, 1, 0.42];
     proportionOfPointsToCrop = 0.03;
%      p = create_manifold_repairing_data_v2(p, peak_point, proportionOfPointsToCrop);
     
     np = size(p,1);
     dim = size(p, 2);  
     %%%%%%%%%%%%%% add noise %%%%%%%%%%%%%%%%
     pp_noiseless = p;
     p1 = ones(size(p, 1), dim);
     noise = 0.1;0.2;%0.1;
   
     p = pp_noiseless + noise*(rand(np, dim) - 0.5*p1); 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q_1 = p(randsample(size(p,1), 200, false), :);
    
%     D = pdist([p; peak_point]);
%     D = squareform(D);
%     selected_pointsDist = D(end,:);
%     [v,s] = sort(selected_pointsDist);
%     c = s(2);%find(i == 1);
%     p_c = p(c,:);    % 6.4817


%     selectedIndexes = s(2:40);
%     q_1 = [q_1; p(selectedIndexes,:)];
    q_1 = unique(q_1,'rows');

    
     % upsample
    if (false)
        nq = size(q_1,1);
        dim = size(q_1, 2);
        p1 = ones(size(q_1, 1), dim);
        qq_1 = q_1;
        num_of_duplicates =1;
        for i=1:num_of_duplicates%10
              noise = 0.2;
            qq_t = q_1 + noise*(rand(nq, dim) - 0.5*p1); % 315x4
            qq_1 = [qq_1;qq_t];
        end
        q_1 = qq_1;
    end
    
    
    q = q_1;

    t = 1;
    niter = 1000;
    DIM_RED = true;
    if (DIM_RED)
        mu = 0;
        sigma = 1;
        new_dim = 20;
        pp_t = p';
        G = normrnd(mu,sigma,[size(pp_t, 2), new_dim]);
        B = pp_t*G;
        [Q,R] = qr(B, 0);
        pp_n = p*Q;

%         if (INITIALISE_DIM_RED_MAT)
            DimRedM = Q;
%         end
    else
        DimRedM = eye(size(p, 2));
    end

    
    figure(1);
    plot3(p(:,1), p(:,2), p(:,3), 'og', q(:,1), q(:,2), q(:,3), 'or'); title('Iter 0');
    view(-110,-20);
    
    t = t+1;   
    R = 9;

    Sp = round(size(p, 1)/size(q, 1));

    %global d;
    d=0.1;

    CALCULATE_H = true;
    if (CALCULATE_H)%CREATE_DATA||DIM_RED)%DIM_RED)%CREATE_DATA)%)
    %     h = 1/3*findHwithMinJPoints(pp, DimRedM, 6);
        minNumOfPoints = 8*Sp;
        
        initial_h_1  = data_obj.getInitialH(p, DimRedM, minNumOfPoints); %8 - elipse;4;%calculateFillDistance(pp, DimRedM, 0.10);*
        

        h1 = calculateFillDistancePerPoint(p, DimRedM, minNumOfPoints, true, initial_h_1, []);

        h1 = median(h1);
        if (sum(isnan(h1)) > 0)
            return;
        end

        index = randsample(size(p,1), round(np/5));
        index = sort(index);
        qq_uniform = p(index,:);
        initial_h_2  = data_obj.getInitialH(qq_uniform, DimRedM, minNumOfPoints/2); 
        
        h2 = calculateFillDistancePerPoint(qq_uniform, DimRedM, minNumOfPoints/2, false, initial_h_2, []);
        h2 = median(h2);
%         h2/h1
%     size(p,1)/size(qq_uniform,1)
%         h2 = h1;
        % h2 = calculateFillDistancePerPoint(pp, DimRedM, minNumOfPoints, false, initial_h, []);

    end

end
%h1 = 0.1;
% h1 = 0.30;h2 = 0.35;

% h_t = 0.15;%3/2;
% h1=h_t; 
% h2 = 0.15;
% h2 = 0.37;
q = q_1;
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_t = 0.10;%3/2;
% debug_point_index = getClosePointIndex(q, [0.5, 0.76, 0.83]) % peak
% debug_point_index = getClosePointIndex(q, [0.78, 0.39, 0.02]) %flat
% debug_point = [0.12, 1, 0.78];
% debug_point_index = getClosePointIndex(q, debug_point);


figure(10);%Ssubplot(3,4,t); 
plot3(p(:,1), p(:,2), p(:,3), 'og', q_1(:,1), q_1(:,2), q_1(:,3), 'or'); 
view(-110,-20); 

view(-110,-20);
kk0 = 1;
c2 = [];
% hole_Radius = 0.45;
% for i=1:size(q,1)
%     peak_weight = norm(q(i,:) - peak_point);
%     c2(i) = theta(peak_weight, hole_Radius);
% end
% 
% c2 = (c2-min(c2))/(max(c2)-min(c2));
% c2=(1-c2);
% 
% 
% figure(1);scatter3( q(:,1), q(:,2),  q(:,3), [],double(c2)+0.01); colorbar;


% h_t_dense = 0.25;
% [c3, minValue, maxValue, hole] = calculateDensityCoef(p, p, DimRedM, d,  h_t_dense, true);
% % figure(1);scatter3( p(:,1), p(:,2),  p(:,3),round(c3*100)+1); colorbar;
% figure(1);scatter3( p(:,1), p(:,2),  p(:,3), double(hole)+0.01); colorbar;


% data_obj.draw_point_by_weight(q, c2 ,  q(1,:), q(1,:), 1, 1);
% view(-110,-20);
% 
% figure(7);close(7);
% draw_point_by_weight(q, c2 + 0.001, q(1,:),  q(1,:), 7, 100, 'c');
% view(-110,-20);


%%%%%%%%%%%%%%%%
for j=1:niter
    
    if (mod(j, 50) == 0)
        j
    end
        
	if (j>1)
        q_k_1 = q_k;
        df_k_1 = df_k;
    end
   
	q_k = q;
     
    for i=1:size(q,1)
        
        odp = repmat(q(i,:), size(p, 1), 1);     
        odq = repmat(q(i,:), size(q, 1), 1);

        eeq = calculateNorm(q(i, :), q, DimRedM);
        eep = calculateNorm(q(i, :), p, DimRedM);
       
%         h1=h_t; h2=h_t;

        %%%%%%%%%%%%%%%% 4. Real derivative %%%%%%%%%  
         Hd_p = Hd(eep,d);
        hdvec = (theta(eep, h1)./Hd_p).*(1 - 2./h1^2 * Hd_p.^2); % the real derive with respect to q,c and Hd

        [eta, d_eta] = eta_function(eeq);

        betvec = 2*(theta(eeq, h2)./eeq).*(abs(d_eta) + 1/h2^2 *eta.*eeq ); % redivative by varios eta, derive with respect to q,
        betvec(i,1)=0;
        
        v = hdvec'*(odp-p)/sum(hdvec);
        u = betvec'*(odq-q);
        %%%%%%%%%%%%%%%% END Real derivative %%%%%%%%%  
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         h2=h_t;
%         betvec2 = 2*(theta(eeq, h2)./eeq).*(abs(d_eta) + 1/h2^2 *eta.*eeq );
%         betvec2(i,1)=0;
%         sparsity = (sum(abs(betvec2)> 0.01));
%         sp(i) = 1;
        
       

%         peak_weight = norm(q(i,:) - peak_point);
% 
%         %peak_weight = 0.36
%         coef = exp(-(peak_weight)^2/.8)*5;
%         sp(i) = coef;
%         h3 = h1;h_t*coef*3/4;
        
%         h1 = h_t*abs(4-coef)/3; % let the points a bit more freedom to move from p
%         Hd_p = Hd(eep,d);
        
%         R = 3;
% %         hdvec = (theta(eep, h1).*(Hd_p-R)./Hd_p).*(2 - 2./h1^2 * Hd_p.*(Hd_p-R)); % the real derive with respect to q,c and Hd
% %         v = hdvec'*(odp-p)/sum(hdvec);
%         Hd_p = Hd(eep,d);
%         hdvec = (theta(eep, h1)./Hd_p).*(1 - 2./h1^2 * Hd_p.^2); % the real derive with respect to q,c and Hd
%         v = hdvec'*(odp-p)/sum(hdvec);
%         
        
%         D = pdist(q);
%         Z = squareform(D);
%         L = theta(Z, h3);
%         L(logical(eye(size(L)))) = 0;
%         L(logical(eye(size(L)))) = -sum(L, 2);
%         x_ = L*q;
%         
        e = theta(eeq, h3);
%         x_i = L*(q-odq);%x_i(i,:)
%         
%         smoth_with_q = 2*sum((q-odq).*e).*sum((2/h3^2*(q-odq).^2-1).*e); %1/peak_weight.^2
        
        smoth_with_q=1;
        coef = 1;

        if (j==1)
            norm_v = sqrt(v*v');
            norm_u = sqrt(u*u');
            norm_smoth_with_q = sqrt(smoth_with_q*smoth_with_q');
            w1(i) = 1/norm_v; %1/sum(hdvec);%
            w2(i) = 2*1/norm_u; %1/sum(betvec);% 
            w3(i) = 3*1/(sum(e.^2));%(2/h3^2-1)%1/norm_smoth_with_q;
           
            a(i) = abs(sum(hdvec));
            w(i) = sign(sum(hdvec))*norm_v/norm_u; %sign(sum(hdvec))* ;*SHIRA
            
        else
        end

        
        %%%%%%%%%%%%%%%% Update the gamma Gradient Descent coeficient    %%%%%%%%%%%%%%%%%%%%%%
 
%         df_k(i,:) = norm_v*(- w2(i)*u + w1(i)*v + w3(i)*smoth_with_q );%
%         df_k(i,:) = v - w2(i)/w1(i)*u;

%         w_q = w2(i)/w1(i);% w2(i)/(a(i)*w1(i));
% %         df_k(i,:) =   (v - w_q*u);%1/ a(i) * + w3(i)*smoth_with_q;%;%+
% 
%         c1 = (peak_weight-0.2)^2;%(((coef-4)^4)/256);
        df_k(i,:) =  w1(i)*v  - w2(i)*u ; %
        
%         df_k(i,:) = 1/(sum(hdvec))*(v - w(i)*u);% 1/(sum(hdvec))* SHIRA


        
        if (j<2)
            gama(i) = 0.1; %0.01;
        else
            d_q = q(i,:) - q_k_1(i,:);
            d_deriv = df_k(i,:) - df_k_1(i,:);
        
            gama(i) = (d_q*d_deriv')/(d_deriv*d_deriv');
            if ( isnan(gama(i)))
                 gama(i) = 0;
                 ['const gamma ' num2str(i) ' ' num2str(j)]
            end
        end
        
        %%%%%%%%%%%%%%%% q^(k+1) = q^(k) + gamma*df/dq(E1+E2) %%%%%%%%%%%
        next_q_i= q(i,:) - gama(i)*df_k(i,:);
        
        nanq = isnan(next_q_i);
        if (sum(nanq(:)) > 0)
           'eeror'
             figure(2);plot3(p(:,1), p(:,2), p(:,3), '*b', q(:,1), q(:,2), q(:,3), '*r', q(i,1), q(i,2), q(i,3), '*g');   
        end
       

        
        q(i,:) = next_q_i;
     
%         if (i==1)
%             figure(2);%subplot(3,4,t); 
%             plot3(p(:,1), p(:,2), p(:,3), 'b', q(:,1), q(:,2), q(:,3), '*r', q(1,1), q(1,2), q(1,3), '*g'); 
%         end
       
        med_data_min_points_q(kk0, 1:7) = [w1(i)*median(v), w2(i)*median(u), w3(i)*median(smoth_with_q), sum(abs(hdvec)>0.01), sum(abs(betvec)> 0.01), sum(abs(e)> 0.01) , coef];
        kk0 = kk0 + 1;
    end
    
    kk0 = 1;
    
    if (j==1)
        w_temp = w;
%         figure(2);plot(w_temp)
        w_shift = 0;
        w_coef = 1;
        
        plot(1:size(w1, 2), w1, 'r', 1:size(w2, 2), w2);%, 'g', 1:size(w3, 2), w3, 'b');
        w1 = w_coef*median(w1)*ones(size(w1))+w_shift;%15*ones(size(w));10
        w2 = w_coef*median(w2)*ones(size(w2))+w_shift;
        w3 = w_coef*median(w3)*ones(size(w3))+w_shift;
        
        a = w_coef*median(a)*ones(size(a));
%         w4 = w_coef*median(w4)*ones(size(w4))+w_shift;
        
    end

    med_data_min_points(kk1, 1:6) = median(med_data_min_points_q(:, 1:6));%[w1(i)*median(v), w2(i)*median(u),  sum(abs(hdvec)>0.01), sum(abs(betvec)> 0.01), sum(abs(L(:,1))> 0.01)];
    kk1 = kk1+1;
    
%   figure(2); plot(p(:,1), p(:,2), '*b', q(:,1), q(:,2), '*r');

   
    
    if (mod(j, 10) == 0)
%        [ + w1(i)*v; - w2(i)*u; w3(i)*smoth_with_q]
       figure(2); plot3(p(:,1), p(:,2), p(:,3), 'og', q(:,1), q(:,2), q(:,3), 'or');
       view(-110,-20); 
%         [hd, relative_err_mean, relative_err_max] = data_obj.evaluateError(p, q, DimRedM);
%         figure(3); scatter3(q(:,1), q(:,2), q(:,3), 10, med_data_min_points_q(:,4)); colorbar
%         
%         figure(4); scatter3(q(:,1), q(:,2), q(:,3), 10, sp); colorbar
%     showCurrentPointMove
        'l';
    end

end

'l';