%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function calculate the norm of the difference of the two given
% vectors (p1, and p2), while using linear sketching to avoid the effect of
% noise in high dimension. The linease sketching is performed using the DimRedM matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n_dist = calculateNorm(p1, p2, DimRedM)

odq = repmat(p1, size(p2, 1), 1);

% Perform linear sketching, and reduce the dimension of the vectors using
% the DimRedM matrix
odq_dr = odq * DimRedM;
p2_dr = p2 * DimRedM;

% Calculate the norm in the low dimension
ep = p2_dr-odq_dr;

n_dist = sqrt(sum(abs(ep).^2,2));

