classdef Data_Bunny < Data_Interface
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
            ptCloud = pcread('D:\Google Drive\Phd\LOP\lopmatlabcode_nD\MainCode\examples\3dScaning\bunny\bunny\data\bun000.ply');
            
            pp = ptCloud.Location;
            
            pp = pp(randsample(size(pp,1), 1000, false), :);
            pp = unique(pp,'rows');
       end
              
       
       function lt = evaluateError(obj, pp, qq, DimRedM, dim, innerDim)
            for i=1:size(qq, 1)
                A = qq(i,1:innerDim^2);
                A = reshape(A, innerDim, innerDim);
                M = (A*A' - eye(innerDim));
                lt(i) = norm(M(:));

            end
            
       end
       
       function drawData(obj, pp, qq, figNum, innerDim, c1)
           
           figure(figNum); close(figNum);
           figure(figNum); 
%            pcshow(qq);
           hold on
           scatter3(pp(:,1), pp(:,2), pp(:,3), '.g')             
           scatter3(qq(:,1), qq(:,2), qq(:,3), '.r')             
           hold off
          
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