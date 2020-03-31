%compute asymmetricity of the droplet, difference between droplet tip and top

%close all
asymm = zeros(loop,1);

%loop on time
for i = 1:loop
    
    %disp(i)
    %ite = round(time/deltaT);
    
    K = KKK(:,i);
    %asymm(i) = abs((abs(K(end))-abs(K(1)))/(abs(K(end))+abs(K(1))));
    asymm(i) = abs((K(end)-K(1))/(K(end)+K(1)));
    %asymm(i) = abs(K(end)-K(1));
    
    if i==15000
        break
    end
    
end

figure
plot(0:deltaT:(i-1)*deltaT,asymm(1:i),'r','LineWidth',2)
xlabel('t')
ylabel('|k_{tip}-k_{top}|')
grid on