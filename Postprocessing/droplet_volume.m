%compute droplet volume after every oteration

V = zeros(PARAM.loop,1);

for i = 1:PARAM.loop
    
    disp(i)
    a = risa(PARAM.m+2:end,i+1);
    b = risb(PARAM.m+2:end,i+1);
    
    V(i) = axis_int_gauss_vect(a',b');
    
end

figure
plot(abs(V-V(1))/V(1),'o-')
xlabel('ite')
ylabel('err Volume')
grid on
    
    