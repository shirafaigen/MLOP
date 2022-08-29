%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function calculate the "norm" that is robust to noise
% The input is a vector which we need to calculate its Hd norm, the d is
% the \epsilon parameter for the calculateion
% it returns a vectore of the norm values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hd=Hd(t,d)
hd = sqrt(t.^2 + d^2);
