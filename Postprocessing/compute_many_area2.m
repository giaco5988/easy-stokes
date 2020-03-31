%compute areas at different iterations

I = 1500;
step = 5;
dt = deltaT;

area1 = zeros(numel(1:step:I),1);
area2 = zeros(numel(1:step:I),1);
area3 = zeros(numel(1:step:I),1);
area4 = zeros(numel(1:step:I),1);

% area1(1) = surf_int_spline(risa(1:nbrel(1)+1,1)',risb(1:nbrel(1)+1,1)');
% area2(1) = surf_gauss(risa(1:nbrel(1)+1,1)',risb(1:nbrel(1)+1,1)');
% area3(1) = surf_lailai(risa(1:nbrel(1)+1,1)',risb(1:nbrel(1)+1,1)');
% area4(1) = surf_gauss_cylind(risa(1:nbrel(1)+1,1)',risb(1:nbrel(1)+1,1)');

for i=1:I/step
    
    area1(i) = surf_int_spline(risa(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)');
    area2(i) = surf_gauss(risa(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)');
    area3(i) = surf_lailai(risa(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)');
    area4(i) = surf_gauss_cylind(risa(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel(1+(i-1)*step)+1,1+(i-1)*step)');
    
%take in account the bug
%     area1(i) = surf_int_spline(risa(1:nbrel((i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel((i-1)*step)+1,1+(i-1)*step)');
%     area2(i) = surf_gauss(risa(1:nbrel((i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel((i-1)*step)+1,1+(i-1)*step)');
%     area3(i) = surf_lailai(risa(1:nbrel((i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel((i-1)*step)+1,1+(i-1)*step)');
%     area4(i) = surf_gauss_cylind(risa(1:nbrel((i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel((i-1)*step)+1,1+(i-1)*step)');
    
end

figure
plot(dt:dt*step:I*dt,area1)
hold on
plot(dt:dt*step:I*dt,area2,'or')
plot(dt:dt*step:I*dt,area3,'k')
plot(dt:dt*step:I*dt,area4,'go')
xlabel('t')
ylabel('A')
hold off
legend('cartesian trapezi','cartesian gauss','cylindrical trapezi','cylindrical gauss')