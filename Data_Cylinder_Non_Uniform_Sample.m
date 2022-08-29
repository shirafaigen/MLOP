%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the cylinder structure manifold object
% It contains methods like generate and draw data of the cylinder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Data_Cylinder_Non_Uniform_Sample < Data_Interface
    properties

    end
    
    methods
        
        % Find an initial approximation for the fill-distance parameter,
        % which will later be improved
        function h = getInitialH(obj, p, DimRedM, minNumOfPoints)
            % find the minNumOfPoints closest neighbors, and find their "typical" distance
            N = size(p, 1);

            num_to_sample = round(N*0.20); %calculate this estimate to 20% of the data for fast evaluation
            indexes = (round(rand(round(num_to_sample), 1)*(N-1))+1);%, num_to_sample);
            size(indexes)
            size(p)
            num_to_sample
            for kk=1:num_to_sample
                currentPoint = p(indexes(kk),:);
                eep = calculateNorm(currentPoint, p, DimRedM); 
                [v, pla] = sort(eep);
                h_k(kk) = median(eep(pla(1:minNumOfPoints)));
            end
            h = median(h_k);
        end
        
        function [w_coef, w_shift, h_coef] = getLOPCoefficients(obj)
            w_coef = 1;
            w_shift = 0;
            h_coef = 0;
        end
        

        function folder = getFolder(obj)
            folder = 'Cylinder';
        end


        % Generate the cylinder data which needs to be reconstgructed
        function pp = createData2(obj, dim, innerDim, N, add_noise)
            pp = obj.createData_with_resolution(dim, innerDim, N, add_noise, 0.06); 
        end

        % Generate the cylinder data with specified resolution (step size)
        function pp = createData_with_resolution(obj, dim, innerDim, N, add_noise, step)
            
            % initialize the parametrs
            t = 0:step:2; 
            u = (0.1:step:1.5)*pi; 
            sizet = size(t); 
            sizeu = size(u); 
            np = sizet(2)*sizeu(2);

            p0 = zeros(np, dim);
            pp = p0;
            R=1; 

            
            if(add_noise)
                noise = 0.2;
            else
                noise = 0;
            end
            
            for i=1:sizet(2)
                for j=1:sizeu(2)
                      % calculate the cylinder coords using: pp((i-1)*sizeu(2)+j,:) = t(i)*[1 1 1 1] + R/sqrt(2)*(cos(u(j))*[0 1 -1 0] + sin(u(j))*[1 0 0 -1]);
                      v1 = zeros(1, dim);
                      v1(1:4) = 1;

                      v2 = zeros(1, dim);
                      v2(1:4) = [1 0 0 0];

                      v3 = zeros(1, dim);
                      v3(1:4) = [0 1 0 0];

                      if(add_noise)
                          c_u = u(j) + noise*(rand(1,1) - 0.5*u(j)); 
                          c_t = t(i) + noise*(rand(1, 1) - 0.5*t(i));
                      else
                          c_u = u(j); 
                          c_t = t(i);
                
                      end
                      
                      pp((i-1)*sizeu(2)+j,:) = c_t*v1+ R^2/sqrt(2)*(cos(c_u)*v2 + sin(c_u)*v3);
                      pp((i-1)*sizeu(2)+j,(5:end)) = 0;noise*(rand(size(v1(5:end),2), 1)' - 0.5*v1(5:end));                        
                end
            end
        
       end
              
       function [parameters, metaDataStruct] = getIntrinsicParametersFromData(obj, pp, metaDataStruct)
           parameters = pp(:,1:3);
       end
       
       % Function that draws the weight of each point. Used for debugging
       function draw_point_by_weight(obj,  qq, weight_p, current_q_i, next_q_i, figNum, coef, minima)
                
            figure(figNum); close(figNum);
            figure(figNum)

            weight_p = abs(weight_p)*coef;
            weight_p(weight_p<0.01) = 0.7;

            hold on
            q=qq;
            x=q(:,1);y=q(:,2);z=q(:,3);
            scatter3(x,y,z,20, weight_p);

            scatter3(current_q_i(:,1) , current_q_i(:,2), current_q_i(:,3),'m', 'filled');
            scatter3(next_q_i(:,1) , next_q_i(:,2), next_q_i(:,3),'c', 'filled');

            hold off;
            colorbar
       end
        
       % Plot the cylinder data in 3D. if both P, and Q sets are given with
       % plot both, otherwise just P
       function drawData(obj, pp, qq, figNum, innerDim)
           figure(figNum); close(figNum)
           
            if (isempty(qq))
                figure(figNum);
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'og');
            else
                figure(figNum); hold on
                
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'og');
                
                q=qq;
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'or');
                hold off
            end
           view(-110,-20);
       end
   end
   
end