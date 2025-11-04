% 历时 217.982941 秒。
clear;
clc
close all;

nt=2540*2;
isnap=60;    % snapshot sampling

dx=15;
h=15;
%nx=430;
% nx=450;
% %nz=220;
% nz=450;


dt=0.001; % calculate time step from stability criterion
% dt=0.002; % calculate time step from stability criterion
% dt=0.001; % calculate time step from stability criterion
tau=dt;


%f0=13;

%t0=2/f0;                       % initialize time axis

f0=15*pi;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian

nt=nt-1;
B = ifft(exp(-2i*sin([0:nt-1]*pi/(2*nt))'*[0:nt-1]),2*nt,'symmetric');
T = B(1:nt,1:nt) ; % <- The Forward Time Dispersion Transform matrix
src = T * src(:);    % <- Transforming the 1xN source-time vector f

load('vp.mat')
[nz,nx]=size(dataxx);
vp=dataxx;

load('vs.mat')
vs=dataxx;

nz=nz+100;

vp1=zeros(nz,nx);
vs1=zeros(nz,nx);

for ii=1:200
    for jj=1:nx
        vp1(ii+50,jj)=vp(ii,jj);
    end
end

for ii=1:50  %%top
    for jj=1:nx
        vp1(ii,jj)=vp(1,jj);
    end
end


for ii=nz-49:nz
    for jj=1:nx
        vp1(ii,jj)=vp(end,jj);
    end
end


for ii=1:200
    for jj=1:nx
        vs1(ii+50,jj)=vs(ii,jj);
    end
end

for ii=1:50  %%top
    for jj=1:nx
        vs1(ii,jj)=vs(1,jj);
    end
end


for ii=nz-49:nz
    for jj=1:nx
        vs1(ii,jj)=vs(end,jj);
    end
end

vp=vp1;
vs=vs1;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates

rou=1;

xs=nx/2;

zs=51;

p=zeros([nz nx]); Vx=p; Vz=p;
Txxx=p;
Txzz=p;
Tzzz=p;
Txzx=p;

Txx=p;
Txz=p;
Tzz=p;


Vxx=p;
Vzz=p;
Vxz=p;
Vzx=p;
coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];;

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordTxx=zeros(nt,nx);
seis_recordTxz=zeros(nt,nx);
tic
for it=1:nt-2,
    %Txx/x
    %     Txxx=coeff(1)*(Txx-circshift(Txx,[0 1]))+...
    %         coeff(2)*(circshift(Txx,[0 -1])-circshift(Txx,[0 2]))+...
    %         coeff(3)*(circshift(Txx,[0 -2])-circshift(Txx,[0 3]))+...
    %         coeff(4)*(circshift(Txx,[0 -3])-circshift(Txx,[0 4]))+...
    %         coeff(5)*(circshift(Txx,[0 -4])-circshift(Txx,[0 5]))+...
    %         coeff(6)*(circshift(Txx,[0 -5])-circshift(Txx,[0 6]))+...
    %         coeff(7)*(circshift(Txx,[0 -6])-circshift(Txx,[0 7]));

    Txx1= (Txx)-circshift(Txx,[0,1]);
    Txx2= (circshift(Txx,[0,-1])-circshift(Txx,[0,2]));
    Txx3= ( circshift(Txx,[0,-2])-circshift(Txx,[0,3]));
    Txx4=  ( circshift(Txx,[0,-3])-circshift(Txx,[0,4]));
    Txx5=  ( circshift(Txx,[0,-4])-circshift(Txx,[0,5]));
    Txx6=  ( circshift(Txx,[0,-5])-circshift(Txx,[0,6]));
    Txx7= ( circshift(Txx,[0,-6])-circshift(Txx,[0,7]));


    Txxx=coeff(1)*Txx1+coeff(2)*Txx2+coeff(3)*Txx3+coeff(4)*Txx4+coeff(5)*Txx5+coeff(6)*Txx6+coeff(7)*Txx7;


    %Txz/z
    %     Txzz=coeff(1)*(circshift(Txz,[ 1])-circshift(Txz,[ 0]))+...
    %         coeff(2)*(circshift(Txz,[ 2])-circshift(Txz,[ -1]))+...
    %         coeff(3)*(circshift(Txz,[ 3])-circshift(Txz,[ -2]))+...
    %         coeff(4)*(circshift(Txz,[ 4])-circshift(Txz,[ -3]))+...
    %         coeff(5)*(circshift(Txz,[ 5])-circshift(Txz,[ -4]))+...
    %         coeff(6)*(circshift(Txz,[ 6])-circshift(Txz,[ -5]))+...
    %         coeff(7)*(circshift(Txz,[ 7])-circshift(Txz,[ -6]));
    Txzz1=(circshift(Txz,[ 1])-circshift(Txz,[ 0]));
    Txzz2=(circshift(Txz,[ 2])-circshift(Txz,[ -1]));
    Txzz3=(circshift(Txz,[ 3])-circshift(Txz,[ -2]));
    Txzz4=(circshift(Txz,[ 4])-circshift(Txz,[ -3]));
    Txzz5=(circshift(Txz,[ 5])-circshift(Txz,[ -4]));
    Txzz6=(circshift(Txz,[ 6])-circshift(Txz,[ -5]));
    Txzz7=(circshift(Txz,[ 7])-circshift(Txz,[ -6]));


    Txzz=coeff(1)*Txzz1+coeff(2)*Txzz2+coeff(3)*Txzz3+coeff(4)*Txzz4+coeff(5)*Txzz5+coeff(6)*Txzz6+coeff(7)*Txzz7;

    %Tzz/z
    %     Tzzz=coeff(1)*(circshift(Tzz,[ 0])-circshift(Tzz,[ -1]))+...
    %         coeff(2)*(circshift(Tzz,[ 1])-circshift(Tzz,[ -2]))+...
    %         coeff(3)*(circshift(Tzz,[ 2])-circshift(Tzz,[ -3]))+...
    %         coeff(4)*(circshift(Tzz,[ 3])-circshift(Tzz,[ -4]))+...
    %         coeff(5)*(circshift(Tzz,[ 4])-circshift(Tzz,[ -5]))+...
    %         coeff(6)*(circshift(Tzz,[ 5])-circshift(Tzz,[ -6]))+...
    %         coeff(7)*(circshift(Tzz,[ 6])-circshift(Tzz,[ -7]));

    Tzzz1=(circshift(Tzz,[ 0])-circshift(Tzz,[ -1]));
    Tzzz2=(circshift(Tzz,[ 1])-circshift(Tzz,[ -2]));
    Tzzz3=(circshift(Tzz,[ 2])-circshift(Tzz,[ -3]));
    Tzzz4= (circshift(Tzz,[ 3])-circshift(Tzz,[ -4]));
    Tzzz5=(circshift(Tzz,[ 4])-circshift(Tzz,[ -5]));
    Tzzz6= (circshift(Tzz,[ 5])-circshift(Tzz,[ -6]));
    Tzzz7= (circshift(Tzz,[ 6])-circshift(Tzz,[ -7]));

    Tzzz=coeff(1)*Tzzz1+coeff(2)*Tzzz2+coeff(3)*Tzzz3+coeff(4)*Tzzz4+coeff(5)*Tzzz5+coeff(6)*Tzzz6+coeff(7)*Tzzz7;

    %Txz/x
    %     Txzx=coeff(1)*(circshift(Txz,[0 -1])-Txz)+...
    %         coeff(2)*(circshift(Txz,[0 -2])-circshift(Txz,[0 1]))+...
    %         coeff(3)*(circshift(Txz,[0 -3])-circshift(Txz,[0 2]))+...
    %         coeff(4)*(circshift(Txz,[0 -4])-circshift(Txz,[0 3]))+...
    %         coeff(5)*(circshift(Txz,[0 -5])-circshift(Txz,[0 4]))+...
    %         coeff(6)*(circshift(Txz,[0 -6])-circshift(Txz,[0 5]))+...
    %         coeff(7)*(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));

    Txzx1=(circshift(Txz,[0 -1])-Txz);
    Txzx2=(circshift(Txz,[0 -2])-circshift(Txz,[0 1]));
    Txzx3=(circshift(Txz,[0 -3])-circshift(Txz,[0 2]));
    Txzx4=(circshift(Txz,[0 -4])-circshift(Txz,[0 3]));
    Txzx5=(circshift(Txz,[0 -5])-circshift(Txz,[0 4]));
    Txzx6=(circshift(Txz,[0 -6])-circshift(Txz,[0 5]));
    Txzx7=(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));

    Txzx=coeff(1)*Txzx1+coeff(2)*Txzx2+coeff(3)*Txzx3+coeff(4)*Txzx4+coeff(5)*Txzx5+coeff(6)*Txzx6+coeff(7)*Txzx7;

    Vx=Vx+1./(rou).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rou).*dt.*(Tzzz+Txzx)/h;

    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);


    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,50,50,0.007);

    Vxx=coeff(1)*(circshift(Vx,[0 -1])-Vx)+...
        coeff(2)*(circshift(Vx,[0 -2])-circshift(Vx,[0 1]))+...
        coeff(3)*(circshift(Vx,[0 -3])-circshift(Vx,[0 2]))+...
        coeff(4)*(circshift(Vx,[0 -4])-circshift(Vx,[0 3]))+...
        coeff(5)*(circshift(Vx,[0 -5])-circshift(Vx,[0 4]))+...
        coeff(6)*(circshift(Vx,[0 -6])-circshift(Vx,[0 5]))+...
        coeff(7)*(circshift(Vx,[0 -7])-circshift(Vx,[0 6]));

    Vxz=coeff(1)*(circshift(Vx,[ 0])-circshift(Vx,[ -1]))+...
        coeff(2)*(circshift(Vx,[ 1])-circshift(Vx,[ -2]))+...
        coeff(3)*(circshift(Vx,[ 2])-circshift(Vx,[ -3]))+...
        coeff(4)*(circshift(Vx,[ 3])-circshift(Vx,[ -4]))+...
        coeff(5)*(circshift(Vx,[ 4])-circshift(Vx,[ -5]))+...
        coeff(6)*(circshift(Vx,[ 5])-circshift(Vx,[ -6]))+...
        coeff(7)*(circshift(Vx,[ 6])-circshift(Vx,[ -7]));

    Vzz=coeff(1)*(circshift(Vz,[ 1])-circshift(Vz,[ 0]))+...
        coeff(2)*(circshift(Vz,[ 2])-circshift(Vz,[ -1]))+...
        coeff(3)*(circshift(Vz,[ 3])-circshift(Vz,[ -2]))+...
        coeff(4)*(circshift(Vz,[ 4])-circshift(Vz,[ -3]))+...
        coeff(5)*(circshift(Vz,[ 5])-circshift(Vz,[ -4]))+...
        coeff(6)*(circshift(Vz,[ 6])-circshift(Vz,[ -5]))+...
        coeff(7)*(circshift(Vz,[ 7])-circshift(Vz,[ -6]));

    Vzx=coeff(1)*(Vz-circshift(Vz,[0 1]))+...
        coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
        coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
        coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
        coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
        coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
        coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));

    Txx=Txx+dt*rou.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rou.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rou.*(vs.^2).*(Vxz+Vzx)/h;

    Txx(zs,xs)=Txx(zs,xs)+src(it);
    Tzz(zs,xs)=Tzz(zs,xs)+src(it);

    seis_recordTxx(it,:)=Txx(zs,:);
    seis_recordTxz(it,:)=Txz(zs,:);

    if rem(it,isnap)== 0,
        imagesc(x,z,Vx), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
end
toc
save('figure10f.mat')