function eigentest()
t = zeros(1,20);
s = zeros(1,20);
t_m = zeros(1,20);
for i=1:20
  A = rand(i*5);
  %A = hessenberg(A);
  tic;
  eigenvalues(A);
  t(i) = toc;

  tic;
  eig(A);
  t_m(i) = toc;

  s(i)= i*5;
end

figure
ax1 = subplot(1,2,1);
p1 = plot(s,t,'LineWidth',1);
ax1.XGrid = 'on';
ax1.YGrid = 'on';
p1.Marker = 'o';
title('QR Algorithm Runtime');
ylabel('Time (sec)');
xlabel('Matrix size (nxn)');

ax2 = subplot(1,2,2);
p2 = plot(s,t_m,'r','LineWidth',1);
ax2.XGrid = 'on';
ax2.YGrid = 'on';
p2.Marker = 'o';
title('eig() Runtime');
ylabel('Time (sec)');
xlabel('Matrix size (nxn)');
print('Runtimes without Hessenberg form', '-depsc');

end
