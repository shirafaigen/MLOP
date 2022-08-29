%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function calculate the gaussian function of the given parameter t,
% with sigma equals to h
% it returns a vectore of the gaussian values
% t is a vector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = theta(t, h)
w = exp(-t.^2/h^2);
end
