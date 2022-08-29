%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the interface for the object/manifold which we want to perform the MLOP
% It contains methods like generate and draw data
% Each objec that wants to be reconstructed, need to inherit this interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Data_Interface < handle
   
    properties
       m_refPP % Reference data for error calculation
       m_lowResStep = 0.005; % The step rerolution for ref data creation
       
       m_dim   % The dimension of the data
       
       m_N %Number of points to generate 
    end
   
   % Methods that must be implemented
   methods (Abstract)
       
       % Generate the data of the current manifold
       pp = createData2(obj, dim, innerDim, N, add_noise);
       
       % Extract the intrinsic parameters from the data
       [parameters, metaDataStruct] = getIntrinsicParametersFromData(obj, pp, metaDataStruct);
       
       % Generate data of the current manifold, with certain resolution
       % (step parameter)
       pp = createData_with_resolution(obj, dim, innerDim, N, add_noise, step);

   end
   
    methods
        
        % Default values of the coefficient
        function [w_coef, w_shift, h_coef] = getLOPCoefficients(obj)
            w_coef = 1;
            w_shift = 0;
            h_coef = 3;
        end
        
       % Generate data 
       function pp = createData(obj, dim, innerDim, N, add_noise)
           obj.m_dim = dim;
           obj.m_N = N;
           
           % generate the manifold data, that we would like to reconstruct
           pp = obj.createData2( dim, innerDim, N, add_noise);
           
           % generate the clean, dense manifold data, for error evaluation purposes
           obj.createReferenceData(obj.m_dim, -1, obj.m_N, obj.m_lowResStep);
       end
        
        % Draw the data in parallel coordinates (in case the intrinsic dimension is greater then 3)
        function drawData(obj, pp, qq, figNum, innerDim)
            metaDataStruct = struct('create_bean_bounds', true, 'bean_bounds', []);
           [parameters_p, metaDataStruct] = obj.getIntrinsicParametersFromData(pp, metaDataStruct);
           
           metaDataStruct.create_bean_bounds = false;
           if (~isempty(qq))
            [parameters_q, metaDataStruct] = obj.getIntrinsicParametersFromData(qq, metaDataStruct);
           else
               parameters_q = [];
           end
              
               
           lables_p_q = [repmat({'p'}, size(parameters_p, 1), 1); repmat({'q'}, size(parameters_q, 1), 1)];

           figure(figNum+1);
           parallelcoords([parameters_p; parameters_q], 'Group', lables_p_q);
        end

        % The method creates a reference data on a dense net
        function createReferenceData(obj, dim, innerDim, N, step)
          obj.m_refPP = createData_with_resolution(obj, dim, innerDim, N, false, step);
        end

        % Evaluate the error of the given data vs the reference data m_refPP
        function [hd, relative_err_mean, relative_err_max] = evaluateError(obj, pp, qq, DimRedM, dim, innerDim)
           [hd ] = HausdorffDist(pp,qq);
            
            refPP = obj.m_refPP;
            if (isempty(refPP)) % do not evaluate the error
                err = Inf;
            else
                err = zeros(size(qq,1), 1);
                for i=1:size(qq,1)
                    i
                    eeq = calculateNorm(qq(i, :),  refPP, DimRedM);
                    [v,p] = min(eeq);
                    
                    ref_norm = refPP(p,:) * DimRedM;
                    ref_norm = sqrt(sum(abs(ref_norm).^2,2));

                    err(i) = v/ref_norm;
                end  
               
            end
            relative_err_mean = mean(err);
            relative_err_max = max(err);
            
          [relative_err_mean, relative_err_max]
        end
    end   
end