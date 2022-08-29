classdef Data_long_Cone_Non_Uniform_Sample < Data_Interface
   properties
      m_step
      m_ref_data
   end
   methods
       
       function obj = Data_Cone_Non_Uniform_Sample()
            obj.m_step = 0.08;
       end
         
       function setStep(obj, new_value)
           obj.m_step = new_value;
       end
       
       function setRefData(obj, ref_data)
           obj.m_ref_data = ref_data;
       end
       
%        function pp = createData(obj, dim, innerDim, N, add_noise)
%            pp = createData2(obj, dim, innerDim, N, add_noise, 0.08);
%        end
        function h = getInitialH(obj, pp, DimRedM, minNumOfPoints)
           global ALPHA_BETA_FROM_PAPER;
           if (ALPHA_BETA_FROM_PAPER)
               h = 1.6;
           else
               h=0.5;%0.5;
           end
        end
       
         function [w_coef, w_shift, h_coef] = getLOPCoefficients(obj)
            w_coef = 1;
            w_shift = 0;
%             h_coef = 2;1.5; % dense data
            h_coef = 1.5;
        end
       
       function folder = getFolder(obj)
           folder = 'cone';
       end
%        function pp = createData(obj, dim, innerDim, N, add_noise)
%            pp = createData2(obj, dim, innerDim, N, add_noise, 0.08);
%        end
       
        function pp = createData2(obj, dim, innerDim, N, add_noise)
            pp = obj.createData_with_resolution(dim, innerDim, N, add_noise, 0.09);
        end

        function pp = createData_with_resolution(obj, dim, innerDim, N, add_noise, step)

            step = 0.05;%0.06; 
            
            T=4; 
            t = [0:step:1.40, 1.40:2*step:3, 3:4*step:4.8];  0:step:T; 
            u = (0.1:step:1.5)*pi; 
            sizet = size(t); 
            sizeu = size(u); 
            np = sizet(2)*sizeu(2);


%             p0 = zeros(np, dim);
%             pp = p0;
            maxR = 1.5;
            R = 0:maxR/(size(t, 2)-1):maxR;
            
            if(add_noise)
                noise = 0.2;
            else
                noise = 0;
            end
            
            max_u = round(sizeu(2):-(sizeu(2)-2)/sizet(2):1);
            pp = [];
            for i=1:sizet(2)
                
                for j=1:max_u(i)%1:sizeu(2)%max_u(i)
            %         pp((i-1)*sizeu(2)+j,:) = t(i)*[1 1 1 1] + R/sqrt(2)*(cos(u(j))*[0
            %         1 -1 0] + sin(u(j))*[1 0 0 -1]);
                      v1 = zeros(1, dim);
                      v1(1:4) = 1;

                      v2 = zeros(1, dim);
                      v2(1:4) = [0 1 -1 0];

                      v3 = zeros(1, dim);
                      v3(1:4) = [1 0 0 -1];

                      if(add_noise)
                      c_u = u(j) + noise*(rand(1,1) - 0.5*u(j)); 
                      c_R = R(i) + noise*(rand(1, 1) - 0.5*R(i)); 
                      c_t = t(i) + noise*(rand(1, 1) - 0.5*t(i));
                
                      else
                          c_u = u(j); 
                          c_R = R(i); 
                          c_t = t(i);
                
                      end
                      
%                       pp((i-1)*sizeu(2)+j,:) 
                      p_current = c_t*v1+ exp(-c_R^2)/sqrt(2)*(cos(c_u)*v2 + sin(c_u)*v3);
                      p_current(5:end) = 1 + noise*(rand(size(v1(5:end),2), 1)' - 0.5*v1(5:end));
                      
                      pp = [pp; p_current];
                      
%                       pp((i-1)*sizeu(2)+j,(5:end)) = 1 + noise*(rand(size(v1(5:end),2), 1)' - 0.5*v1(5:end));
                            
%                       pp((i-1)*sizeu(2)+j,:) = t(i)*v1+ R(i)/sqrt(2)*(cos(u(j))*v2 + sin(u(j))*v3);

                             
                end
                       

   
            end
            
%             t = sizet(2):step:size(R,2);
%             
%             for i=sizet(2):size(R,2)
%                 pp((i-1)*sizeu(2)+j,:) = t(i)*v1+ exp(-R(i)^2)/sqrt(2)*(cos(0)*v2 + sin(0)*v3);
%             end
       end
              
       
%        function lt = evaluateError(obj, pp, qq, DimRedM, dim, innerDim)
% %             lt = loptest(qq,dim); 
% %             figure(12);plot(lt,'r');
%             
%             for i=1:size(qq, 1)
%                 if (mod(i, 50) == 0)
%                      i
%                  end
%                 
% %                 for j=1:size(obj.m_ref_data, 1)
%                    currentPoint = qq(i, :);
%                    eep = calculateNorm(currentPoint, obj.m_ref_data, DimRedM).^2; 
%                    [v, p] = min(eep);
%                    lt(i) = v;
% %                 end
%             end
%             
%        end
       
       function [parameters, metaDataStruct] = getIntrinsicParametersFromData(obj, pp, metaDataStruct)
           parameters = pp(:,1:3);
       end
       
        function drawData(obj, pp, qq, figNum, innerDim)
            if (isempty(qq))
                figure(figNum);
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'g');
            else
                figure(figNum)
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'*g', 'MarkerFaceColor','g');
                hold on
                q=qq;
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'*r');
                
%                 k=118;
%                 scatter3(x(k),y(k),z(k),'*r');
% %                 
%                 k=119;
%                 scatter3(x(k),y(k),z(k),'*m');
%                 scatter3(x(35),y(35),z(35),'k');
                hold off
           end
        end
        
        function drawData2(obj, pp, qq, figNum, innerDim, i, minPlace)
                figure(figNum)
                q=pp; 
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'g');
                hold on
                q=qq;
                x=q(:,1);y=q(:,2);z=q(:,3);
                scatter3(x,y,z,'r');
                
                x=q(i,1);y=q(i,2);z=q(i,3);
                scatter3(x,y,z,'b');
                
                
                x=q(minPlace,1);y=q(minPlace,2);z=q(minPlace,3);
                scatter3(x,y,z,'b');
                
                hold off;
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