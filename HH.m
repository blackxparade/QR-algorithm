function A = HH(v)

% param v is a column vector!
%
% Householder matrix formula:
%       2
% I - ----- h*h'
%      h'*h



% STEP 1: make an identity matrix with a correct size.
I = eye(length(v));

% STEP 2: do the math
A = I - ( (2/(v' * v)) * (v * v') );

end
