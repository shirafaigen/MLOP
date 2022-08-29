classdef Data_Orthogonal_Matrix_Uniform_Sample < Data_Interface
   properties
      
   end
   methods
       
       function [parameters, metaDataStruct] = getIntrinsicParametersFromData(obj, pp, metaDataStruct)
           parameters = pp(:,1:3);
       end

    
    
        function h = getInitialH(obj)
           global ALPHA_BETA_FROM_PAPER;
           if (ALPHA_BETA_FROM_PAPER)
               h = 0.4;
           else
               global innerDim;
               if innerDim ==2
                    h = 0.4;
               else
                   h=0.01;0.08;0.01;
               end
           end
        end
       
       function folder = getFolder(obj)
            folder = 'O2';
       end
       
     
       function pp = createData2(obj, dim, innerDim, N, add_noise)
            pp = obj.createData_with_resolution(dim, innerDim, N, add_noise, 0.06);
        end

        function pp = createData_with_resolution(obj, dim, innerDim, N, add_noise, step)
            add_noise = true;
            pp = zeros(N, dim);
            
            if (add_noise)
                 switch innerDim
                     case 2
                        step = (pi-(-pi))/N;
                        thetas = -pi:step:pi;
                        
%                          thetas = rand(N, 1);
%                         thetas = thetas*pi;
%                         sing_teta = randi([0 1], N, 1);
%                         sing_teta(sing_teta==0) = -1;
%                         thetas = thetas.*sing_teta;
                     case 3
                        angle=-pi:0.1:pi;
                        thetas=angle(randi(numel(angle),N, innerDim));

%                         noise = 0.9;
%                         p1 = zeros(N, innerDim);
%                         thetas = thetas + noise*(rand(N, innerDim) - 0.5*p1); 
                 end
            else
                switch innerDim
                     case 2
                        thetas = rand(N, 1);
                        thetas = thetas*pi;
                        sing_teta = randi([0 1], N, 1);
                        sing_teta(sing_teta==0) = -1;
                        thetas = thetas.*sing_teta;
                        
                     case 3
                        angle=-pi:0.1:pi;
                        thetas=angle(randi(numel(angle),N, innerDim));
                end
            end
    
    
            for i=1:N

                switch innerDim
                   %%%%%%%%%%%%%%
                    case 2
%                         theta = randn(1);
                        angle=-pi:0.1:pi;
                        theta=thetas(i); %angle(randi(numel(angle),1,1));
                        
                        
                        % theta = 45 * pi/180; % The rotation angle, converted to radians
                        
                        c = cos(theta);
                        s = sin(theta);
                        
                        A = [c  -s ; s  c];
%                         A = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];

                    case 3
                        
%                         A = vrrotvec2mat(thetas);
                        Rx = rotx(thetas(i, 1));
                        Ry = roty(thetas(i, 2));
                        Rz = rotz(thetas(i, 3));
                        A = Rz*Ry*Rx;

                        
                        
                        
                    otherwise
                        A = orth(randn(innerDim, innerDim));

                end

                pp(i, 1:size(A,1)*size(A,2)) = A(:);
            end
            

        switch innerDim
                   %%%%%%%%%%%%%%
                case 2    
                noise = 0.1;
                
        
                case 3    
                noise = 0.2;0.0001;
        end
        
        np = size(pp,1);
        dim = size(pp, 2);
        p1 = ones(size(pp, 1), dim);
% p1(:,1:9) = 0;
n = dim;10
        pp(:,n:end) = pp(:,n:end) + noise*(rand(np, dim-n+1) - 0.5*p1(:,n:end)); % 315x4

       end
              
       
       function lt = evaluateError(obj, pp, qq, DimRedM, dim, innerDim)
            for i=1:size(qq, 1)
                A = qq(i,1:innerDim^2);
                A = reshape(A, innerDim, innerDim);
                M = (A*A' - eye(innerDim));
                lt(i) = norm(M(:));

            end
            
       end
       
       function draw_point_by_weight(obj, p, weight_p, current_p1, current_p2, fignum, coef, title_c)

            figure(fignum);hold on;
            hold on
            scatter( p((weight_p>=0),1), p((weight_p>=0),2), abs(weight_p(weight_p>=0))*coef,'r')
            
            scatter( p((weight_p<0),1), p((weight_p<0),2), abs(weight_p(weight_p<0))*coef,'db')

            scatter( current_p1(:,1), current_p1(:,2), 50,'om', 'filled')
            scatter( current_p2(:,1), current_p2(:,2), 50,'oc','filled')
            
            
            hold off;
       end
       function drawData(obj, pp, qq, figNum, innerDim, c1)
           
           figure(figNum); close(figNum);
           if (isempty(qq))
               values = pp;
           else
               values = qq;
           end
           innerDim = 2;
           switch innerDim
                %%%%%%%%%%%%%%
                case 2
                %%%%%%%%%%%%%%
                    
%                     for i=1:size(values, 1)
%                         A = values(i,1:innerDim^2);
%                         x(i) = A(1,1);
%                         y(i) = A(1,2);
% 
%                     end
%                     
%                     if (isempty(qq))
%                         figure(figNum); viscircles([0 0], 1, 'EdgeColor','g'); hold on;
%                         plot(x, y, 'x'); hold off
%                     else
%                         
%                         for i=1:size(pp, 1)
%                             A = pp(i,1:innerDim^2);
%                             x_p(i) = A(1,1);
%                             y_p(i) = A(1,2);
%                         end
%                     
                        
                        figure(figNum);
                        viscircles([0 0], 1, 'EdgeColor','b'); hold on;
                        
%                         scatter(x, y , 20, c1);
                    
%                         plot(x_p, y_p, '*g',  x,y, '*r'); hold off %, 'markersize',20
                        %plot(x_p, y_p, '*g', x,y, '*r'); hold off
                      plot(pp(:, 1), pp(:, 2), '*g', qq(:, 1), qq(:, 2), '*r'); hold off
                        
                        
%                     end
                    
                %%%%%%%%%%%%%%
                case 3
                %%%%%%%%%%%%%%
                    [x, y, z] = obj.getPointsFromAngleR3(values, innerDim);
                   
                    
                    if (isempty(qq))
                         
                        figure(figNum); [x_s,y_s,z_s] = sphere;
                        scatter3(x_s(:),y_s(:),z_s(:),'g'); hold on
                        scatter3(x,y,z,'r'); hold off
                    else
                        

                        figure(figNum); 
                        [x_s,y_s,z_s] = sphere;
%                         a = gradient(z_s)*0+0.2;
                        hold on
                        scatter3(x_s(:),y_s(:),z_s(:), '.b')
%                         surf(x_s,y_s,z_s,'AlphaData',a, 'FaceAlpha','flat',  'FaceColor','red'); 
                        
                        [x_s, y_s, z_s] = obj.getPointsFromAngleR3(pp, innerDim);
                      
                        scatter3(x_s(:),y_s(:),z_s(:),'.g');
                        
                        scatter3(x,y,z,'r'); 
                        scatter3(x(1),y(1),z(1),'*r'); 

%                        figure(figNum);  scatter3(x_s(:),y_s(:),z_s(:),'g'); hold on
%                         scatter3(x,y,z,'r'); hold off
                        
                    end
                    
           end
          
       end
       
       function [x, y, z] = getPointsFromAngleR3(obj, values, innerDim)
           
            for i=1:size(values, 1)
                A = values(i,1:innerDim^2);
                A = reshape(A, innerDim, innerDim);

                % https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
                % t = acos(0.5*(A(1,1)+A(2,2) + A(3,3)-1));
                % x(i) = (A(3,2)+A(2,3))/(2*sin(t));
                % y(i) = (A(1,3)+A(3,1))/(2*sin(t));
                % z(i) = (A(2,1)+A(1,2))/(2*sin(t));

                r = vrrotmat2vec(A);
                x(i) = real(r(1));
                y(i) = real(r(2));
                z(i) = real(r(3));                     

            end
       end
                    
   end
   
end
%                         y(i) = asin(A(1,3));
%                         z(i) = acos(A(1,1)/cos(y(i)));
%                         x(i) = acos(A(3,3)/cos(y(i)));
%  pitch = asin(A(1,3));
%                         yaw = acos(A(1,1)/cos(pitch));
%                         roll  = asin(A(3,2)/cos(pitch));
% x(i) = -cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll);
%                         y(i) = -sin(yaw)*sin(pitch)*sin(roll)+cos(yaw)*cos(roll);
%                         z(i) =  cos(pitch)*sin(roll);


% x(i) = cos(yaw)*cos(pitch);
% y(i) = sin(yaw)*cos(pitch);
% z(i) = sin(pitch);

                    %   sy = A(1,1) * A(1,1) +  A(2,1) * A(2,1);
% 
%                         singular = sy < 1e-6;
% 
%                         if  (~singular)
%                             x(i) = atan2(A(3,2) , A(3,3));
%                             y(i) = atan2(-A(3,1), sy);
%                             z(i) = atan2(A(2,1), A(1,1));
%                         else 
%                             x(i) = atan2(-A(2,3), A(2,2));
%                             y(i) = atan2(-A(3,1), sy);
%                             z(i) = 0;
%                         end