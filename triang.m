function [Q,R] = triang(A)

[m,n] = size(A);
Q = eye(n);

for i=1:n-1
	B = eye(n);
	B(i:end, i:end) = HH( norm(A(i:end, i), 2) * eye(n+1-i,1) - A(i:end, i));
	Q = Q * B;
	A = B * A;
	 
end
R = A;
end
