function EIG = eigenvalues(A)

[m,n] = size(A);
A = hessenberg(A);

for i=1:10000
	[Q,R] = triang(A);
	A = R * Q;
end

disp(A);

end
