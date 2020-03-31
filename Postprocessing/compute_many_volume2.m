%compute areas at different iterations

I = 41;
%step = 10;
dt = deltaT;

V1 = zeros(numel(I),1);
V2 = zeros(numel(I),1);

% V1(1) = axis_int(risa(1:nbrel(1)+1,1)',risb(1:nbrel(1)+1,1)');
% V2(1) = axis_int_gauss(risa(1:nbrel(1)+1,1)',risb(1:nbrel(1)+1,1)');

for i = 1:I
    
    V1(i) = axis_int(Y(i,1:2:end-1),Y(i,2:2:end));
    V2(i) = axis_int_gauss(Y(i,1:2:end-1),Y(i,2:2:end));
    
%take in account the bug
%     V1(i) = axis_int(risa(1:nbrel((i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel((i-1)*step)+1,1+(i-1)*step)');
%     V2(i) = axis_int_gauss(risa(1:nbrel((i-1)*step)+1,1+(i-1)*step)',risb(1:nbrel((i-1)*step)+1,1+(i-1)*step)');
    
end

figure
plot(T,V1/4/pi*3)
hold on
plot(T,V2/4/pi*3,'or');
xlabel('t')
ylabel('Volume')
hold off
legend('cartesian trapezi','cartesian gauss')