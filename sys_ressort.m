close all
clear all;

dt = 1e-8;	% Pas d’intégration en s
Niter=1e6;	% Nombre d’itérations
k0 = 2e-8;	% Raideur des ressorts en N/m
uma = 1.66e-27;	% Unité de masse atomique en kg
m = 204*uma ;	% Masse des polymere en kg
vtrac = 1e-2;	% Vitesse traction en m/ s
theta=0.5;
Natom=7;
K=k0/m*[2 -1 0 0 0 0 ; -1 2 0 0 -1 0 ; 0 0 2 -1 0 0 ; 0 0 -1 2 -1 0 ; 0 -1 0 -1 3 -1 ; 0 0 0 0 -1 2 ];
B=[zeros(6), K; -eye(6), zeros(6)];
Y=[0.02;0.02;0;0;0.01;0.01;0.01]
V=zeros(12,Niter+1);
D = inv(eye(12)+dt*theta*B);
for i=1:Niter
    V(:,i+1)=D*( (eye(12)-dt*(1-theta)*B)*V(:,i)+[0;0;0;0;0;k0*vtrac*dt*(i-0.5)/m;0;0;0;0;0;0]);
    X=[V(7:12,i+1)'*dt,vtrac*dt*(i+1)];
    if(mod(i,1000)==0)
        plot(X,Y,'.r','MarkerSize',10);
        axis([0 vtrac*dt*Niter -0.1 0.1]);
        hold;
        drawnow;
    end
end