function HESS = hessenberg(A)

[m,n] = size(A);
HESS = A;

for i=1:n-2
	column = norm(HESS(i+1:end,i),2) * eye(length(HESS(i+1:end,i)),1) - HESS(i+1:end,i);
	HH = eye(length(column)) - 2/(column' * column) * (column * column');
	B = eye(n);
	%disp(size(HH));
	%disp(size(B));
	B(i+1:end, i+1:end) = HH;
	HESS = B * HESS * B;
end
end
