%when computing quantities in the domain with BEM, move point really close
%to the interface on the closest node of the interface and assign the
%corrisponding quantity of the node

function [X0,Y0,U,V,p,InOut] = close_interface_visu(a,b,X0,Y0,U,V,p,InOut,solution,crit)

    count = 0;
    
    crit_dist = max(sqrt(diff(a).^2+diff(b).^2))/crit;
    for i = 1:numel(X0)
            
            dist = sqrt((X0(i)-a).^2+(Y0(i)-b).^2);
            
            if min(dist) < crit_dist
                
                [~,k] = min(dist);
                
                X0(i) = a(k);
                Y0(i) = b(k);
                U(i) = solution(2*k-1);
                V(i) = solution(2*k);
                %p(i) = nan;
                %p(i) = (p(i-1)+p(i+1))/2;
                %ppp = (x())
                
%                 figure(1)
%                 hold on
%                 plot(Y0(i),-X0(i),'or')
%                 quiver(Y0(i),-X0(i),V(i),-U(i),'k')
%                 hold off
                
                InOut(i) = 1;
                
                %check usage
                count = count+1;
                
            end
            
    end
    
    display(['replace interrogation point ' num2str(count) ' times'])
    
end