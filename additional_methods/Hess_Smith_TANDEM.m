clear all
close all
clc

%% Select bodies
    % M  : ordinata massima della linea media             (una cifra)
    % P  : posizione sulla corda dell'ordinata massima M  (una cifra)
    % SS : spessore massimo del profilo simmetrico        (due cifre)
    %
    % n  : numero di punti lungo la corda distribuiti uniformemente,
    %      inclusi i punti del bordo di attacco e del bordo di uscita,  
    %      in corrispondenza dei quali si calcolano i punti sul dorso
    %      e sul ventre
    % M, P e SS possono essere forniti come NUMERI INTERI, compresi, 
    %      rispettivamente,
    %
    %      M  fra 1 e 9,  in  %  della corda      
    %      P  fra 1 e 9   in 1/10  della corda   
    %      SS fra 1 e 99  in  %  della corda
    %
    %      oppure, se si preferisce, come VALORI REALI  < 1  
    %      frazione della corda
    %
    % xd and yd are the coordinate of the top. "d" stands for "dorso"
    % xv and yv are the coordinate of the bototm. "v" stands for "ventre"
% 1
M1 = 0; P1=0; SS1=10; n1=100; % NACA0010
%2
M2 = 0; P2=0; SS2=8; n2=100; % NACA0008
%
x12 = 1.2;  % distance between forward airfoil LE and back airfoil LE, chord length normalized to 1
y12 = 0;

%% WIND DATA 
U = 1; % asymptotic wind
alfa1 = 2*pi/180; alfa2 = -2*pi/180; % angle of attack of the airfoils

Uinf(1)=U*cos(alfa1); Uinf(2)=U*sin(alfa1);

%% LOAD POINTS: Panel extrema
 [xd1, yd1, xv1, yv1] = NACA_4d (M1, P1, SS1, n1, 'halfcos');
 [xd2, yd2, xv2, yv2] = NACA_4d (M2, P2, SS2, n2, 'halfcos');
 
 figure(901)
 plot(xd1,yd1,'r',xv1,yv1,'r')
 hold on
 plot(xd2,yd2,'b',xv2,yv2,'b')
 axis equal
 title('Generated airfoils shape')
 
 xe1=[flipud(xv1(:,1)); xd1(2:end,1)];   
 ye1=[flipud(yv1(:,1)); yd1(2:end,1)]; 
 Np1=size(xe1,1)-1;
 
 xie2 =[flipud(xv2(:,1)); xd2(2:end,1)];   
 etae2=[flipud(yv2(:,1)); yd2(2:end,1)]; 
 Np2=size(xie2,1)-1;
 
 gamma = alfa2 - alfa1;
 
 xie2 =   xie2 - 0.5*(xd2(end,1)+xd2(1,1));
 xe2  =   xie2 * cos(gamma) + etae2 * sin(gamma); 
 ye2  = - xie2 * sin(gamma) + etae2 * cos(gamma); 
 xe2  =   xe2  + 0.5*(xd2(end,1)+xd2(1,1));
 
 xe2 = xe2 + x12; 
 ye2= ye2 + y12;

 figure(902)
 plot(xe1,ye1,'r+')
 hold on
 plot(xe2,ye2,'b+')
 axis equal
 title('Positioning of airfoils in space')

%% Build geometry information needed
% Control points, each point is in the middle of the panel
xp1 = 0.5*(xe1(1:Np1) + xe1(2:Np1+1)); 
yp1 = 0.5*(ye1(1:Np1) + ye1(2:Np1+1));

xp2 = 0.5*(xe2(1:Np2) + xe2(2:Np2+1)); 
yp2 = 0.5*(ye2(1:Np2) + ye2(2:Np2+1));

figure(902)
plot(xp1,yp1,'ko')
hold on
plot(xp2,yp2,'ko')
axis equal
title('Positioning of airfoils in space with control points')
legend('Airfoil 1','Airfoil 2','Control points 1','Control points 2')

% normals
normal1 = zeros(Np1,2); normal2 = zeros(Np2,2); % In column 1 we put x component, in column 2 the y component
normal1(:,1) = -(ye1(2:Np1+1)-ye1(1:Np1))./(sqrt((ye1(2:Np1+1)-ye1(1:Np1)).^2+(xe1(2:Np1+1)-xe1(1:Np1)).^2));
normal1(:,2)=   (xe1(2:Np1+1)-xe1(1:Np1))./(sqrt((ye1(2:Np1+1)-ye1(1:Np1)).^2+(xe1(2:Np1+1)-xe1(1:Np1)).^2)); 
normal2(:,1) = -(ye2(2:Np2+1)-ye2(1:Np2))./(sqrt((ye2(2:Np2+1)-ye2(1:Np2)).^2+(xe2(2:Np2+1)-xe2(1:Np2)).^2));
normal2(:,2)=   (xe2(2:Np2+1)-xe2(1:Np2))./(sqrt((ye2(2:Np2+1)-ye2(1:Np2)).^2+(xe2(2:Np2+1)-xe2(1:Np2)).^2)); 
% tangent vector
tangent1 = zeros(Np1,2); tangent2 = zeros(Np2,2);
tangent1(:,1) = (xe1(2:Np1+1)-xe1(1:Np1))./(sqrt((ye1(2:Np1+1)-ye1(1:Np1)).^2+(xe1(2:Np1+1)-xe1(1:Np1)).^2));
tangent1(:,2) = (ye1(2:Np1+1)-ye1(1:Np1))./(sqrt((ye1(2:Np1+1)-ye1(1:Np1)).^2+(xe1(2:Np1+1)-xe1(1:Np1)).^2));
tangent2(:,1) = (xe2(2:Np2+1)-xe2(1:Np2))./(sqrt((ye2(2:Np2+1)-ye2(1:Np2)).^2+(xe2(2:Np2+1)-xe2(1:Np2)).^2));
tangent2(:,2) = (ye2(2:Np2+1)-ye2(1:Np2))./(sqrt((ye2(2:Np2+1)-ye2(1:Np2)).^2+(xe2(2:Np2+1)-xe2(1:Np2)).^2));
% Lengths
length1= sqrt((ye1(2:Np1+1)-ye1(1:Np1)).^2+(xe1(2:Np1+1)-xe1(1:Np1)).^2);
length2= sqrt((ye2(2:Np2+1)-ye2(1:Np2)).^2+(xe2(2:Np2+1)-xe2(1:Np2)).^2);

figure(903)
plot(dot(normal1,tangent1),'r')
title('Scalar product between normals and tangents')

%% Build the linear system Ax=b
A11 = zeros(Np1,Np1); A22 = zeros(Np2,Np2);
A12 = zeros(Np1,Np2); A21 = zeros(Np2,Np1);

a11 = zeros(Np1,1); a22 = zeros(Np2,1);
a12 = zeros(Np1,1); a21 = zeros(Np2,1);

C11 = zeros(1,Np1); C22 = zeros(1,Np2);
C12 = zeros(1,Np2); C21 = zeros(1,Np1);

c11 = 0; c22=0; 
c12 = 0; c21=0;

b=zeros(Np1+Np2+2,1);

% build A11 and a11
for i = 1: Np1
    for j = 1: Np1
        [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(i),yp1(i)); % panel in the global reference frame
        [vor.u,vor.v] = vortex( xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(i),yp1(i));  
        A11(i,j) = sou.u*normal1(i,1) + sou.v*normal1(i,2);
        a11(i,1) = a11(i,1) + vor.u*normal1(i,1) + vor.v*normal1(i,2);
    end
end
% build A22 and a11
for i = 1: Np2
    for j = 1: Np2
        [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(i),yp2(i));
        [vor.u,vor.v] = vortex( xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(i),yp2(i));  
        A22(i,j) = sou.u*normal2(i,1) + sou.v*normal2(i,2);
        a22(i,1) = a22(i,1) + vor.u*normal2(i,1) + vor.v*normal2(i,2);
    end
end
% build A12 and a12
for i = 1: Np1
    for j = 1: Np2
        [sou.u,sou.v] = source( xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(i),yp1(i));
        [vor.u,vor.v] = vortex( xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(i),yp1(i));  
        A12(i,j) = sou.u*normal1(i,1) + sou.v*normal1(i,2);
        a12(i,1) = a12(i,1) + vor.u*normal1(i,1) + vor.v*normal1(i,2);
    end
end
% build A21 and a21
for i = 1: Np2
    for j = 1: Np1
        [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(i),yp2(i));
        [vor.u,vor.v] = vortex( xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(i),yp2(i));  
        A21(i,j) = sou.u*normal2(i,1) + sou.v*normal2(i,2);
        a21(i,1) = a21(i,1) + vor.u*normal2(i,1) + vor.v*normal2(i,2);
    end
end
% build C11 and c11
for j = 1:Np1
    [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(1),yp1(1)); 
    [vor.u,vor.v] = vortex(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(1),yp1(1));
    C11(1,j) = sou.u*normal1(1,2)-sou.v*normal1(1,1);
    c11 = c11 + vor.u*normal1(1,2) - vor.v*normal1(1,1);
    [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(Np1),yp1(Np1));    
    [vor.u,vor.v] = vortex(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(Np1),yp1(Np1)); 
    C11(1,j) = C11(1,j) +  sou.u*normal1(Np1,2)-sou.v*normal1(Np1,1);
    c11 = c11 + vor.u*normal1(Np1,2) - vor.v*normal1(Np1,1);    
end
% build C22 and c22
for j = 1:Np2
    [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(1),yp2(1)); 
    [vor.u,vor.v] = vortex(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(1),yp2(1));
    C22(1,j) = sou.u*normal2(1,2)-sou.v*normal2(1,1);
    c22 = c22 + vor.u*normal2(1,2) - vor.v*normal2(1,1);
    [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(Np2),yp2(Np2));    
    [vor.u,vor.v] = vortex(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(Np2),yp2(Np2)); 
    C22(1,j) = C22(1,j) +  sou.u*normal2(Np2,2)-sou.v*normal2(Np2,1);
    c22 = c22 + vor.u*normal2(Np2,2) - vor.v*normal2(Np2,1);    
end
% build C12 and c12
for j = 1:Np2
    [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(1),yp1(1)); 
    [vor.u,vor.v] = vortex(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(1),yp1(1));
    C12(1,j) = sou.u*normal1(1,2)-sou.v*normal1(1,1);
    c12 = c12 + vor.u*normal1(1,2) - vor.v*normal1(1,1);
    [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(Np1),yp1(Np1));    
    [vor.u,vor.v] = vortex(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(Np1),yp1(Np1)); 
    C12(1,j) = C12(1,j) +  sou.u*normal1(Np1,2)-sou.v*normal1(Np1,1);
    c12 = c12 + vor.u*normal1(Np1,2) - vor.v*normal1(Np1,1);    
end
% build C21 and c21
for j = 1:Np1
    [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(1),yp2(1)); 
    [vor.u,vor.v] = vortex(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(1),yp2(1));
    C21(1,j) = sou.u*normal2(1,2)-sou.v*normal2(1,1);
    c21 = c21 + vor.u*normal2(1,2) - vor.v*normal2(1,1);
    [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(Np2),yp2(Np2));    
    [vor.u,vor.v] = vortex(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(Np2),yp2(Np2)); 
    C21(1,j) = C21(1,j) +  sou.u*normal2(Np2,2)-sou.v*normal2(Np2,1);
    c21 = c21 + vor.u*normal2(Np2,2) - vor.v*normal2(Np2,1);    
end

A=[A11 A12 a11 a12;
      A21 A22 a21 a22;
      C11 C12 c11 c12;
      C21 C22 c21 c22];
  
  
for i = 1: Np1
    b(i) = -Uinf(1)*normal1(i,1) - Uinf(2)*normal1(i,2);
end
for i = 1: Np2
    b(Np1+i) = -Uinf(1)*normal2(i,1) - Uinf(2)*normal2(i,2);
end
b(Np1+Np2+1) = -Uinf(1)*tangent1(1,  1)   - Uinf(2)*tangent1(1,  2) ...
                   -Uinf(1)*tangent1(Np1,1) - Uinf(2)*tangent1(Np1,2);
b(Np1+Np2+2) =  -Uinf(1)*tangent2(1,  1)   - Uinf(2)*tangent2(1,  2) ...
                            -Uinf(1)*tangent2(Np2,1) - Uinf(2)*tangent2(Np2,2);


%% Solve the linear system
sol = linsolve(A,b);

sigma1 = sol(1:Np1);
sigma2 = sol(Np1+1:Np1+Np2);
gamma1 = sol(Np1+Np2+1);
gamma2 = sol(Np1+Np2+2);

%% Compute velocity on control points and pressure ditribution
velu1= zeros(Np1,1); velv1 = velu1;
velu2= zeros(Np2,1); velv2 = velu2;
% airfoil 1
for i =1:Np1
    velu1(i) = Uinf(1); velv1(i)=Uinf(2);
    for j = 1:Np1
        [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(i),yp1(i));
        [vor.u,vor.v] = vortex(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp1(i),yp1(i));
       velu1(i) = velu1(i) + sigma1(j)*sou.u + gamma1*vor.u;
       velv1(i) = velv1(i) + sigma1(j)*sou.v + gamma1*vor.v;
    end
    for j = 1:Np2
        [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(i),yp1(i));
        [vor.u,vor.v] = vortex(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp1(i),yp1(i));
       velu1(i) = velu1(i) + sigma2(j)*sou.u + gamma2*vor.u;
       velv1(i) = velv1(i) + sigma2(j)*sou.v + gamma2*vor.v;
    end
end
% arifoil 2
for i =1:Np2
    velu2(i) = Uinf(1); velv2(i)=Uinf(2);
    for j = 1:Np1
        [sou.u,sou.v] = source(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(i),yp2(i));
        [vor.u,vor.v] = vortex(xe1(j),ye1(j),xe1(j+1),ye1(j+1),xp2(i),yp2(i));
       velu2(i) = velu2(i) + sigma1(j)*sou.u + gamma1*vor.u;
       velv2(i) = velv2(i) + sigma1(j)*sou.v + gamma1*vor.v;
    end
    for j = 1:Np2
        [sou.u,sou.v] = source(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(i),yp2(i));
        [vor.u,vor.v] = vortex(xe2(j),ye2(j),xe2(j+1),ye2(j+1),xp2(i),yp2(i));
       velu2(i) = velu2(i) + sigma2(j)*sou.u + gamma2*vor.u;
       velv2(i) = velv2(i) + sigma2(j)*sou.v + gamma2*vor.v;
    end
end

if (max(velu1.*normal1(:,1) + velv1.*normal1(:,2))>10^(-14)) || (max(velu2.*normal2(:,1) + velv2.*normal2(:,2))>10^(-14))
    disp('There is a bug in the program!')
    stop
end

Vt1 = velu1.*tangent1(:,1) + velv1.*tangent1(:,2);
Vt2 = velu2.*tangent2(:,1) + velv2.*tangent2(:,2);

Cp1 = 1-Vt1.^2/(Uinf(1)^2+Uinf(2)^2);
Cp2 = 1-Vt2.^2/(Uinf(1)^2+Uinf(2)^2);

%%
% Cl = sum (Cp_i * L_i* cos(theta_i) )
 Cl1 = -Cp1'*( length1.*normal1(:,2) );
disp('Cl2'); Cl2 = -Cp2'*( length2.*normal2(:,2) );
disp('Cl1:'); disp(Cl1);
disp('Cl2:'); disp(Cl2);

 figure(2); hold on; box on; grid on;
 plot(xp1(1:Np1/2),-Cp1(1:Np1/2),'b--','LineWidth',2);
 plot(xp1(Np1/2+1:end),-Cp1(Np1/2+1:end),'b:','LineWidth',2);
 plot(xp2(1:Np2/2)-xp2(Np2/2),-Cp2(1:Np2/2),'r--','LineWidth',2);
 plot(xp2(Np2/2+1:end)-xp2(Np2/2),-Cp2(Np2/2+1:end),'r:','LineWidth',2);
 legend box off;
 tx=xlabel('$x$');
 ty=ylabel('$-C_p$');
 tz = legend('bottom airfoil 1','top airfoil 1','bottom airfoil 2','top airfoil 2');
 tx.Interpreter='latex';
 ty.Interpreter='latex';
 tz.Interpreter='latex';
 set(gca,'TickLabelInterpreter', 'latex');
 x0=10;
 y0=10;
 tx.Interpreter='latex';
 ty.Interpreter='latex';
 tz1.Interpreter='latex';
 tz2.Interpreter='latex';
 set(gca,'TickLabelInterpreter', 'latex');
 tz.Interpreter='latex';
 width=600;
 height=400;
 set(gcf,'position',[x0,y0,width,height])
 ax = gca;
 ax.FontSize = 18;

 

%%

% alfa = [ -5 -2.5 0 2.5 5];
% xfoil_CL = [-0.6602 -0.3304 0 0.3304 0.6602];
% CL = [-0.6555 -0.3290 0 0.3290 0.6555];
% 
% figure(); hold on; box on; grid on
% plot(alfa,CL,'kd-','MarkerFaceColor','g')
% plot(alfa,xfoil_CL,'k*')
% legend('present','xfoil')
% legend boxoff
% tx=xlabel('$C_{\ell}$')
% ty=ylabel('$\alpha$')
% tx.Interpreter='latex';
% ty.Interpreter='latex';
% tz.Interpreter='latex';
% set(gca,'TickLabelInterpreter', 'latex');
% x0=10;
% y0=10;
% tx.Interpreter='latex';
% ty.Interpreter='latex';
% tz1.Interpreter='latex';
% tz2.Interpreter='latex';
% set(gca,'TickLabelInterpreter', 'latex');
% tz.Interpreter='latex';
% width=600;
% height=400;
% set(gcf,'position',[x0,y0,width,height])
% ax = gca;
% ax.FontSize = 18;


%% plot velocity field
Nfx=100; Nfy = 100;
x=linspace(-0.5,3,Nfx);
y=linspace(-0.5,0.5,Nfy);

Ufield= zeros(Nfx,Nfy); Vfield = Ufield; body = Ufield;

ydplus1 = interp1(xd1(:,1),yd1(:,1),x);
yvplus1 = interp1(xv1(:,1),yv1(:,1),x);
ydplus2 = interp1(xp2(Np2/2+1:end),yp2(Np2/2+1:end),x);
yvplus2 = interp1(flipud(xp2(1:Np2/2)),flipud(yp2(1:Np2/2)),x);

for i = 1:Nfx
    for j = 1:Nfy
        if (x(i)<xd1(end,1) && x(i)>xd1(1,1) && y(j)<ydplus1(i) && y(j)>yvplus1(i)) || ...
           (x(i)<max(xe2) && x(i)>min(xe2) && y(j)<ydplus2(i) && y(j)>yvplus2(i))
           body(i,j)=1;
        else
            body(i,j)=0;
        end
    end
end
            
for i =1:Nfx
    for j = 1:Nfy
        if body(i,j)==0
            Ufield(i,j) = Uinf(1); Vfield(i,j) = Uinf(2);
           for ip = 1:Np1
              [sou.u,sou.v] = source(xe1(ip),ye1(ip),xe1(ip+1),ye1(ip+1),x(i),y(j));
              [vor.u,vor.v] = vortex(xe1(ip),ye1(ip),xe1(ip+1),ye1(ip+1),x(i),y(j));
              Ufield(i,j) = Ufield(i,j) + sigma1(ip)*sou.u + gamma1*vor.u;
              Vfield(i,j) = Vfield(i,j) + sigma1(ip)*sou.v + gamma1*vor.v;
           end
           for ip = 1:Np2
              [sou.u,sou.v] = source(xe2(ip),ye2(ip),xe2(ip+1),ye2(ip+1),x(i),y(j));
              [vor.u,vor.v] = vortex(xe2(ip),ye2(ip),xe2(ip+1),ye2(ip+1),x(i),y(j));
              Ufield(i,j) = Ufield(i,j) + sigma2(ip)*sou.u + gamma2*vor.u;
              Vfield(i,j) = Vfield(i,j) + sigma2(ip)*sou.v + gamma2*vor.v;
           end
        end            
    end
end
% %%
[X,Y] = meshgrid(x,y);
% figure(3); hold on; box on;
% contourf(X',Y',Ufield,100,'LineStyle','None');
% colormap(flipud(hot));
% plot(xe1,ye1,'Linewidth',2,'Color','k');
% plot(xe2,ye2,'Linewidth',2,'Color','k');
% colorbar('off');
% view(2);
% tx=xlabel('$x$');
% ty=ylabel('$y$');
% tx.Interpreter='latex';
% ty.Interpreter='latex';
% set(gca,'TickLabelInterpreter', 'latex');  
% x0=10;
% y0=10;
% width=1250; 
% height=580;
% set(gcf,'position',[x0,y0,width,height]);
% ax = gca;
% axis equal;
% axis off;
% ax.FontSize = 40;
% title('U Velocity field','Interpreter','latex')
% %%
% figure(4); hold on; box on;
% contourf(X',Y',Vfield,100,'LineStyle','None');
% colormap(bluewhitered(256)); colorbar;
% c=colorbar('northoutside');
% c.TickLabelInterpreter='latex';
% colorbar('off');
% plot(xe1,ye1,'Linewidth',2,'Color','k');
% plot(xe2,ye2,'Linewidth',2,'Color','k');
% view(2);
% tx=xlabel('$x$');
% ty=ylabel('$y$');
% tx.Interpreter='latex';
% ty.Interpreter='latex';
% set(gca,'TickLabelInterpreter', 'latex');  
% x0=10;
% y0=10;
% width=1250; 
% height=580;
% set(gcf,'position',[x0,y0,width,height]);
% ax = gca;
% axis equal;
% axis off;
% ax.FontSize = 40;
% title('V Velocity field','Interpreter','latex')

%%
figure(5); hold on; box on;
Velocity=sqrt(Vfield.^2+Ufield.^2);
contourf(X',Y',Vfield,100,'LineStyle','None');
colormap(bluewhitered(256)); colorbar;
c=colorbar('northoutside');
c.TickLabelInterpreter='latex';
colorbar('off');
plot(xe1,ye1,'Linewidth',2,'Color','k');
plot(xe2,ye2,'Linewidth',2,'Color','k');
view(2);
tx=xlabel('$x$');
ty=ylabel('$y$');
tx.Interpreter='latex';
ty.Interpreter='latex';
set(gca,'TickLabelInterpreter', 'latex');  
x0=10;
y0=10;
width=1250; 
height=580;
set(gcf,'position',[x0,y0,width,height]);
ax = gca;
axis equal;
axis off;
ax.FontSize = 40;
title('Velocity field','Interpreter','latex')
%%
figure(6); hold on; box on;
quiver(X',Y',Ufield,Vfield);
colormap(bluewhitered(256)); colorbar;
c=colorbar('northoutside');
c.TickLabelInterpreter='latex';
colorbar('off');
plot(xe1,ye1,'Linewidth',2,'Color','k');
plot(xe2,ye2,'Linewidth',2,'Color','k');
view(2);
tx=xlabel('$x$');
ty=ylabel('$y$');
tx.Interpreter='latex';
ty.Interpreter='latex';
set(gca,'TickLabelInterpreter', 'latex');  
x0=10;
y0=10;
width=1250; 
height=580;
set(gcf,'position',[x0,y0,width,height]);
ax = gca;
axis equal;
axis off;
ax.FontSize = 40;
title('Velocity field','Interpreter','latex')
%% Source and vortex functions
function [u,v] = source( xe1, ye1, xe2, ye2, x, y) % It seems ok
l=sqrt((xe2-xe1)^2+(ye2-ye1)^2);
xl2 = l; yl2=0;
xl=(x-xe1)*(xe2-xe1)/l+(y-ye1)*(ye2-ye1)/l; yl=-(x-xe1)*(ye2-ye1)/l+(y-ye1)*(xe2-xe1)/l;
r1=sqrt(xl^2+yl^2); r2=sqrt((xl-xl2)^2+(yl-yl2)^2);
theta1=atan2(yl,xl); theta2=atan2(yl-yl2,xl-xl2);
if (abs(theta1)<10^(-12) && abs(theta2)>3); theta1=0; theta2=pi; end
if (abs(theta2)<10^(-12) && abs(theta1)>3); theta2=0; theta1=-pi; end
vl=1/(2*pi)*(theta2-theta1); ul=1/(2*pi)*log(r1/r2);
u=ul*(xe2-xe1)/l-vl*(ye2-ye1)/l; v=ul*(ye2-ye1)/l+vl*(xe2-xe1)/l;
end

function [u,v] = vortex( xe1, ye1, xe2, ye2, x, y) % it seems ok
l=sqrt((xe2-xe1)^2+(ye2-ye1)^2);
xl2 = l; yl2=0;
xl=(x-xe1)*(xe2-xe1)/l+(y-ye1)*(ye2-ye1)/l; yl=-(x-xe1)*(ye2-ye1)/l+(y-ye1)*(xe2-xe1)/l;
r1=sqrt(xl^2+yl^2); r2=sqrt((xl-xl2)^2+(yl-yl2)^2);
theta1=atan2(yl,xl); theta2=atan2(yl-yl2,xl-xl2);
if (abs(theta1)<10^(-12) && abs(theta2)>3); theta1=0; theta2=pi; end
if (abs(theta2)<10^(-12) && abs(theta1)>3); theta2=0; theta1=-pi; end
ul=1/(2*pi)*(theta2-theta1); vl=-1/(2*pi)*log(r1/r2);
u=ul*(xe2-xe1)/l-vl*(ye2-ye1)/l; v=ul*(ye2-ye1)/l+vl*(xe2-xe1)/l;
end
