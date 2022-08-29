%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function calculate the eta function, and its derivative of a given vector
% Using the formula: \eta = 1/(3r^2)
% The input is a vector which we need to calculate its eta, and d eta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta, d_eta] = eta_function(eeq)

%  Calculate eta
eta = 1./(3*eeq.^3);

% Calculate its derivatice
d_eta = -1/(eeq.^4);
d_eta = d_eta';

end