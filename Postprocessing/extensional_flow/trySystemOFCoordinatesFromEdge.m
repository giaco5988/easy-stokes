%verify which is the good system of coordinates able to describe a certain shape

clear variables
close all

%choose which system
polar = 1;
fourier = 0;

%physical parameters
Ca = 0.04;
lambda = 1;

if Ca==0.1
pathUpload = '~/Documents/C++/edge_state_extensional/resultsServer2/';
elseif Ca==0.04
pathUpload = '~/Documents/C++/edge_state_extensional/resultsCaSmall/';
end

%nodes
n = 200;

%legendre modes
modes = 100;
symmetric = 1;
symmetricP = 2;
symmetricPP = 1;

%derivatives unifirm spacing
D = finiteDifference1D(n+1,[2 1],2);
D1 = D(:,:,1);  D2 = D(:,:,2);

%parameters
IDshape = 5;
IDdelta = 16;
Round = 1;

%directory where to upload the results
if numel(num2str(Ca))==4
        UploadResults = [pathUpload 'Ca=' num2str(Ca) '0000_lambda=' num2str(lambda) '.000000_IDshape=' num2str(IDshape) '_IDdelta=' num2str(IDdelta) '_Round=' num2str(Round) '_elem=' num2str(n) '_dt=0.010000_loop=25000_RK=2_CPUs=1'];
    elseif numel(num2str(Ca))==3
        UploadResults = [pathUpload 'Ca=' num2str(Ca) '00000_lambda=' num2str(lambda) '.000000_IDshape=' num2str(IDshape) '_IDdelta=' num2str(IDdelta) '_Round=' num2str(Round) '_elem=' num2str(n) '_dt=0.010000_loop=25000_RK=2_CPUs=1'];
end
name = [UploadResults '/drop0.txt'];
    
try
    A = importdata(name);
catch
    warning('No data')
end

aIN = A.data(:,1)'; bIN = A.data(:,2)';
    
%place in the center
xcm = center_mass(aIN,bIN);
aIN = aIN-xcm;
a = aIN;    b = bIN;

figure
plot([aIN flip(aIN)],[bIN -flip(bIN)])
grid on
axis equal
xlabel('x')
ylabel('r')

%independent varaibles
theta = atan(b./a)';
theta = theta + pi*(theta<0);
 
%try polar coordinates
if polar==1
    
    %metrics term
    dThetadT = D1*theta;
    
    %radius
    r = sqrt(a.^2+b.^2)';
    %r = r(2:end-1);
    %r = ones(numel(theta),1);
    
    %first derivative
    rp = (D1*r)./dThetadT;
    
    %second derivative
    rpp = (D2*r)./dThetadT.^2;
    
    figure
    plot(theta,r,'-')
    hold on
    plot(theta,rp,'-')
    plot(theta,rpp,'-')
    grid on
    xlabel('\theta')
    ylabel('r,rp,rpp')
    legend('r','rp','rpp','Location','Best')
    title('polar coordinates')
    
    %fit legendre polynomia to the data
    f = LegendreSeriePolar(theta,r,modes,symmetric);
    fp = LegendreSeriePolar(theta,rp,modes,symmetricP);
    fpp = LegendreSeriePolar(theta,rpp,modes,symmetricPP);

    figure
    loglog(abs(f),'o-')
    grid on
    xlabel('n')
    ylabel('c_n')
    title('modes coeff r')
    
    hold on
    loglog(abs(fp),'o-')
    grid on
    xlabel('n')
    ylabel('c_n')
    title('modes coeff rp')
    
    
    loglog(abs(fpp),'o-')
    grid on
    xlabel('n')
    ylabel('c_n')
    title('modes coeff rpp')

end

if fourier==1
    
    x = a;  y = b;
    
    figure
    plot(theta,x,'-')
    hold on
    plot(theta,y,'-')
    grid on
    xlabel('\theta')
    ylabel('r,rp,rpp')
    legend('x(\theta)','r(\theta)','Location','Best')
    title('cartesian coordinates')

    %fourier coeff
    %bn = fourierSin(theta,a',modes);
    %an = fourierCos(theta,a',modes);
    
    bn = fft(a(1:end-1));
    %an = fft(a(1:end-1));
    an = fft(a);

    %THETA = [theta; theta+pi]; THETA = THETA(1:end-1);
    %AAA = [a'; flip(a)'];   AAA = AAA(1:end-1);
    %an = fourierCos(THETA,AAA,modes);

%     figure
%     plot([theta; theta+pi],[a'; flip(a)'])
%     grid on

    figure
    %loglog(abs(an(2:2:end)))
    loglog(real(an))
    
    hold on
    %loglog(abs(bn))
    grid on
    xlabel('n')
    ylabel('a_n')
    title('shape')
    
    figure
    loglog(abs(bn(1:2:end)))
    hold on
    %loglog(abs(bn))
    grid on
    xlabel('n')
    ylabel('b_n')
    title('shape')

end













