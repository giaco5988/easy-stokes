%compute series of legendre form given data points

function [f,SumP] =ChebSerieSpectral(theta,x,modes,symmetric,weight,der)

    % differentiation and integration
    %WG = SPECTRAL.WG; d
    
    %initialize variables
    f = zeros(1,numel(modes));
    
    %compute modes corfficiens
    SumP = 0;
    for i = 1:modes

            %I take only symmetric modes if the shape is symmetric
            temp = i-1;
            if symmetric==1
                polN = 2*temp;
            elseif symmetric==0
                polN = temp;
            elseif symmetric==2
                polN = 2*temp+1;
            end
            
%             if der==1
%                 polN = polN+4;
%             elseif der==2
%                 polN = polN+4;
%             end

            %compute coeff
            int = x.*cos(polN*theta);
            f(i) = 2/pi*weight'*int;

            if polN==0
                f(i) = f(i);
            end
            
            %compute serie
            SumP = SumP + f(i)*cos(polN*theta);
            
    end
    
%     figure
%     plot(theta,r)
%     hold on
%     plot(theta,SumP,'--')
%     grid on
%     hold off
    
    %Vserie = 2/3*INT*(r.^3.*sin(theta));
    %err = abs(Vserie-V0)/V0;

end