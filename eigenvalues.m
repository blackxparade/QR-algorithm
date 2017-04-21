function EIG = eigenvalues(A)

[m,n] = size(A);
EIG = eye(n,1);

for i=1:500
	[Q,R] = triang(A);
	A = R * Q;
end

for i=1:n
	EIG(i) = A(i,i);
end

%disp(A);

end
