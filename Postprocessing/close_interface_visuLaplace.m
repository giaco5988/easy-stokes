%when computing quantities in the domain with BEM, move point really close
%to the interface on the closest node of the interface and assign the
%corrisponding quantity of the node

function [X0,Y0,PHI,InOut] = close_interface_visuLaplace(x,y,X0,Y0,PHI,InOut,solution,maxElem,crit)

    count = 0;
    
    crit_dist = maxElem*crit;
    for i = 1:numel(X0)
            
            dist = sqrt((X0(i)-x).^2+(Y0(i)-y).^2);
            
            if min(dist) < crit_dist
                
                [~,k] = min(dist);
                
                X0(i) = x(k);
                if Y0(i)~=0% don't move this point radially if it is located on the axis
                    Y0(i) = y(k);
                end
                PHI(i) = solution(k);
                
                InOut(i) = 1;
                
                %check usage
                count = count+1;
                
            end
            
    end
    
    display(['replace interrogation point ' num2str(count) ' times'])
    
end