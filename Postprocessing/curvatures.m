%test mehod to compute the volume and the surface of a solid, both
%axisymmetric

close all
%clear all

%color for the graphs
color = ['b'; 'r'; 'k'; 'm'; 'g'];
%color = ['-k'; '--k'; 'o-k'];
color = repmat(color,100,1);

%deform = [5 8 11];


check = 1:10:2000;   
k1_south = zeros(numel(check),1);
k1_east = zeros(numel(check),1);
k2_south = zeros(numel(check),1);
k2_east = zeros(numel(check),1);
k_south = zeros(numel(check),1);
k_east = zeros(numel(check),1);

for i = 1:numel(check)

          disp(i)

          %compute the spline coeff
          [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (risa(:,check(i))', risb(:,check(i))');

          %compute the versor normal tp the node (N) and to the element (r)
          %[N,r] =normal_versor_clean(a,b,q-1);
          N = [by./sqrt(bx.*bx+by.*by) by(end)+2*cy(end)+3*dy(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
              -bx./sqrt(bx.*bx+by.*by) -bx(end)-2*cx(end)-3*dx(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

          %%%%%%%%%%%%%%%%%%% COMPUTE CURVATURE WITH SPLINES %%%%%%%%%%%%%%%%%%%%%%

          K1 = curv_spline2(bx,by,cx,cy,dx,dy);

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          %those because the points are on the axis
          N(:,1) = [1; 0];
          N(:,end) = [-1; 0];

          %second component of the curvature
          K2 = N(2,:)./risb(:,check(i))';
          %this because on the axis the I cannot use the definition of before
          K2(1) = K1(1);
          K2(end) = K1(end);

          %sum the two component of the curvature
          K = K1+K2;
          
%           figure(1)
%           plot(K2)

          k1_south(i) = K1(round(numel(K1)/2));%min(K1);
          k1_east(i) = K1(1);%max(K1);
          
          k2_south(i) = K2(round(numel(K1)/2));%min(K2);
          k2_east(i) = K2(1);%max(K2);
          
          k_south(i) = K(round(numel(K)/2));
          k_east(i) = K(1);     

end

figure(1)
hold on
plot(k1_south,k1_south+k1_east,'Linewidth',2)
%plot(k1_east,k1_south+k1_east,'Linewidth',2)
hold off
xlabel('k1_s')
ylabel('k1_s+k1_e')        
grid on
   
figure(2)
hold on
plot(k2_south,k2_south+k2_east,'Linewidth',2)
hold off
xlabel('k2_s')
ylabel('k2_s+k2_e')        
grid on

figure(3)
hold on
plot(k_south,k_south+k_east,'Linewidth',2)
hold off
xlabel('k_s')
ylabel('k_s+k_e')        
grid on


% for k = 1:numel(deform)
%     
%     filename = strcat('Relaxation_q=200_visc=1_Dt=0.1_loop=3000_DELTA=',num2str(deform(k)),'_Ca=10_RK2.mat');
%     load(filename)
% 
%     check = 1:10:2000;
%     k_south = zeros(numel(check),1);
%     k_east = zeros(numel(check),1);
% 
%     for i = 1:numel(check)
% 
%           disp(i)
% 
%           %compute the spline coeff
%           [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (risa(:,check(i))', risb(:,check(i))');
% 
%           %compute the versor normal tp the node (N) and to the element (r)
%           %[N,r] =normal_versor_clean(a,b,q-1);
%           N = [by./sqrt(bx.*bx+by.*by) by(end)+2*cy(end)+3*dy(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
%               -bx./sqrt(bx.*bx+by.*by) -bx(end)-2*cx(end)-3*dx(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];
% 
%           %%%%%%%%%%%%%%%%%%% COMPUTE CURVATURE WITH SPLINES %%%%%%%%%%%%%%%%%%%%%%
% 
%           K1 = curv_spline2(bx,by,cx,cy,dx,dy);
% 
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%           %those because the points are on the axis
%           N(:,1) = [1; 0];
%           N(:,end) = [-1; 0];
% 
%           %second component of the curvature
%           K2 = N(2,:)./risb(:,check(i))';
%           %this because on the axis the I cannot use the definition of before
%           K2(1) = K1(1);
%           K2(end) = K1(end);
% 
%           %sum the two component of the curvature
%           K = K1+K2;
% 
%           k1_south(i) = K1(round(numel(K1)/2));
%           k1_east(i) = K1(1);
%           
%           k2_south(i) = K2(round(numel(K1)/2));
%           k2_east(i) = K2(1);
%           
%           k_south(i) = K(round(numel(K1)/2));
%           k_east(i) = K(1); 
% 
%     end
%     
%     figure(1)
%     hold on
%     plot(k1_south,k_south+k1_east,color(k),'Linewidth',2)
%     hold off
%     
%     figure(2)
%     hold on
%     plot(k2_south,k2_south+k2_east,color(k),'Linewidth',2)
%     hold off
%         
%     figure(3)
%     hold on
%     plot(k_south,k_south+k_east,color(k),'Linewidth',2)
%     hold off
% 
% end
% 
% figure(1)
% xlabel('k1_s')
% ylabel('k1_s+k1_e')        
% grid on
% legend(strcat('D=',num2str(deform(1))),strcat('D=',num2str(deform(2))),strcat('D=',num2str(deform(3))))
% 
% figure(2)
% xlabel('k2_s')
% ylabel('k2_s+k2_e')        
% grid on
% legend(strcat('D=',num2str(deform(1))),strcat('D=',num2str(deform(2))),strcat('D=',num2str(deform(3))))
% 
% figure(3)
% xlabel('k_s')
% ylabel('k_s+k_e')        
% grid on
% legend(strcat('D=',num2str(deform(1))),strcat('D=',num2str(deform(2))),strcat('D=',num2str(deform(3))))
    
%add the ideal path thorugh many ellipses
% ks = k1_south(1):0.01:1;
% 
% figure(3)
% hold on
% plot(ks,1./ks+ks,'--k','Linewidth',2)
% hold off