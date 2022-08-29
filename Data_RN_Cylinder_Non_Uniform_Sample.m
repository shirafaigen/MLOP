classdef Data_RN_Cylinder_Non_Uniform_Sample < Data_Interface
   properties
     m_t
   end
   
    methods
        function h = getInitialH(obj)
           h = 0.4;
        end

        function d = getD(obj)
            d = 20;
        end


        function [w_coef, w_shift, h_coef] = getLOPCoefficients(obj)
            
            d = obj.getD();
            
            switch d
                case 2
                    w_coef = 1;
                    w_shift = 0;
                    h_coef = 1.5;
                    
                case 3
                    w_coef = 0;
                    w_shift = 10;
                    h_coef = 2;
                    
                case 10
                    w_coef = 0;
                    w_shift = 10;
                    h_coef = 3;
                case 20
                    w_coef = 30;
                    w_shift = 0;
                    h_coef = 2;
                    
                    
                otherwise
                    w_coef = 1;
                    w_shift = 0;
                    h_coef = 3;
            end
                
          

        end

        function folder = getFolder(obj)
        folder = 'Cylinder-Rn';
        end
       
        function  t = getT(obj)
           t =  obj.m_t;
        end
        
        function pp = createData2(obj, dim, innerDim, N, add_noise)
          obj.m_t = 0:0.06:2; 
          pp = obj.createData_with_resolution(dim, innerDim, N, add_noise, 0.06);
        end

  

        function pp = createData_with_resolution(obj, dim, innerDim, N, add_noise, step)

            d = getD(obj);

            t = 0:step:2;

            u = (0.1:step:1)*pi; 
            U = repmat(u, d, 1);


            sizet = size(t); 
            sizeu = size(U); 
            np = sizet(2)*sizeu(2);

            pp = zeros(np, dim);
            R=1; 


            if(add_noise)
                noise = 0.2;
            else
                noise = 0; 
            end


            % prepare orthogonal bases
            I = eye(d+1);

            for i=1:sizet(2)
                for j=1:sizeu(2)

                    %%%%%%%%%%%%%%%%% calculating: p_i = t*[1,1,1,1] + (cos+sin)*mul(sin)*v_i %%%%%%%%%%%%%%                  

                    %%%%%%%%%%%%%%%%% First Term :t*[1,1,1,1] %%%%%%%%%%%%%%%%%
                    v = zeros(1, dim);
                    v(1:d+1) = 1;
                    c_t = t(i) + noise*(rand(1, 1) - 0.5*t(i));
                    pp((i-1)*sizeu(2)+j,:) = c_t*v;


                    %%%%%%%%%%%%%%%%% Second Term: (cos+sin)*mul(sin)*v_i %%%%%%%%%%%%%%%%%
                    % Calculations based on: https://en.wikipedia.org/wiki/N-sphere
                    % x1 = R cos(u1)
                    % x2 = R sin(u1) cos(u2)
                    % x3 = R sin(u1) sin(u2) cos(u3)
                    % ...
                    % x_n-1 = R sin(u1) sin(u2) ... sin(u_n-1) cos(u_n)
                    % xn = R sin(u1) sin(u2) ... sin(u_n-1) sin(u_n)

                    sinCosSum = zeros(1, dim);

                    % for every teta_i there should be a net
                    for k=1:d-1

                        sinProd = 1;
                        for m = 1:k-1
                            teta = U(m,j) + noise*(rand(1,1) - 0.5*U(m, j));
                            sinProd = sinProd*sin(teta);
                        end 

                        teta = U(k,j) + noise*(rand(1,1) - 0.5*U(k, j));

                        v = zeros(1, dim);
                        v(1:size(I(:, k), 1)) = I(:, k);

                        % the last element is always cos(u), except for the n'th iteration
                        sinProd_c = cos(teta)*sinProd*v;
                        sinCosSum = sinCosSum + sinProd_c;
                    end 

                    v = zeros(1, dim);
                    v(1:size(I(:, k+1), 1)) = I(:, k+1);

                    % add a sin multiplication
                    sinProd_s = sin(teta)*sinProd*v;
                    sinCosSum = sinCosSum + sinProd_s;

                    % verification: sum(sinCosSum.^2) == 1  (with no noise)
                    pp((i-1)*sizeu(2)+j,:) = pp((i-1)*sizeu(2)+j,:) + R^2*sinCosSum; 

                    % add noise to the other coords
                    pp((i-1)*sizeu(2)+j,(d+2:end)) = 1 + noise*(rand(size(v(d+2:end),2), 1)' - 0.5*v(d+2:end));                         
                end
            end
        end

        % Extract the intrinsic parameters from the data
        function  [parameters, metaDataStruct] = getIntrinsicParametersFromData(obj, pp, metaDataStruct)

            [parameters, metaDataStruct.bean_bounds] = obj.extractAngles(pp, ~metaDataStruct.create_bean_bounds, metaDataStruct.bean_bounds);
        end

        % Extract the teta_i, and t_i from the given vector representation
        % using the following observation: p_i = t*[1,1,1,1] + (cos+sin)*mul(sin)*v_i
        function [tetas, bead_bounds] = extractAngles(obj, pp, is_qq, bead_bounds)
           d = getD(obj);

           %%%%%%%%%%%%%%%%% extract the T's based on the following quadratic equation: -nt^2 + 2sum(p_i)t+(1-sum(p_i^2) = 0
           n = d+1;
           sum_x = sum(pp(:,1:n), 2);
           sum_x_2 = sum(pp(:,1:n).^2, 2);

           b = 2*sum_x;
           a = -n;
           c = 1-sum_x_2;

           t_1 = (-b + sqrt(b.^2-4.*a.*c))./(2*a);
           t_2 = (-b - sqrt(b.^2-4.*a.*c))./(2*a);

           t_1 = (real(t_1));
           t_2 = (real(t_2));

           % since there is noise, the t's are not exactly what they
           % should be. Therefore, devide the range in to size(t2,2) beans,
           % and but every t, to the right bean

           t2 = obj.getT();

           % use only the positive beans, the negavives will be put to the
           % first bean (just a descition)
           pos_index = t_1>0;
           pos_t = t_1(pos_index);


           if (~is_qq)
                % the hist function opens a window, we close it afterwards
                figure(100);
                [b] = histogram(pos_t, size(t2,2));

                bead_bounds = b.BinEdges;
                close(100);
           end

           % for each t, loacate the bean it belong to
           pos_t_bean_num = discretize(pos_t, bead_bounds);

           % incase there is a pot_t> max bean_bound, just use the last bean for this t
           pos_t_bean_num(isnan(pos_t_bean_num)) = size(bead_bounds, 2);

           t(pos_index) = bead_bounds(pos_t_bean_num);
           t(~pos_index) = bead_bounds(1);
           
           % finally, the t's were located extract them from p, to get the
           % sin*cos multiplication
           ts = repmat(t, n, 1)';
           pp(:,1:n) = pp(:,1:n) - ts;
           
           % use the following observation to calculate the tetas one after
           % another. from the first one till the end.
           %        teta = acos(p_i./sin_Multiplication (from 1:i-1))
           for i=1:d-1

               %sin mult
               sinMult = 1;
               for m = 1:i-1
                   sinMult = sinMult.*sin(tetas(:,m));
               end

               tetas(:,i) = real(acos((pp(:,i))./sinMult));
           end
           % since its possiable for teta, to be bigger than pi, normalize
           % it between [0, pi]
           tetas = mod(tetas, pi);

           % concatinate the found t
           tetas = [tetas, t'];
        end



        % Draw the Data, if in d==2 as a cylinder, and also in parallel coords
        function drawData(obj, pp, qq, figNum, innerDim)

           d = getD(obj);

           if (d == 2)
            if (isempty(qq))
                figure(figNum);
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'g');
            else
                figure(figNum)
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'g');
                hold on
                q=qq;
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'b');

                hold off
            end
           end

          % draw in parallel coords
          drawData@Data_Interface(obj, pp, qq, figNum, innerDim)

        end

        function draw_point_by_weight(obj, pp, qq, weight_p, current_q_i, next_q_i, figNum, coef, minima)
                figure(figNum); close(figNum);
                figure(figNum)
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
        %                 scatter3(x,y,z, 'g');
                weight_p = abs(weight_p)*coef;

                weight_p(weight_p<0.01) = 0.7;

                scatter3(x,y,z, weight_p,'g');
                hold on
                q=qq;
                x=q(:,1);y=q(:,2);z=q(:,3);
        %                 scatter3(x,y,z,'r');

                x=qq(current_q_i, 1);y=qq(current_q_i, 2);z=qq(current_q_i, 3);
                scatter3(x,y,z,'b');


        %                 x=q(minPlace,1);y=q(minPlace,2);z=q(minPlace,3);
        %                 scatter3(x,y,z,'b');

                hold off;
        end
        
   end
   
end