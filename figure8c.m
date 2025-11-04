% 历时 41.794250 秒。

clear;
clc
close all;

nt=1000;
isnap=20;    % snapshot sampling

dx=20;
h=20;
%nx=430;
nx=450;
%nz=220;
nz=450;

x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0015; % calculate time step from stability criterion
tau=dt;


%f0=13;
t=(1:nt)*dt;
%t0=2/f0;                       % initialize time axis

f0=16*pi;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


for i=1:nz
    for j=1:nx
        %vp(i,j)=4000;
        vp(i,j)=2500*sqrt(3);
    end
end


for i=1:nz
    for j=1:nx
        %vs(i,j)=2310;
        vs(i,j)=2500;
    end
end

vp(1:nz/2+60,:)=2100*sqrt(3);

vs(1:nz/2+60,:)=2100;
rou=1;

xs=nx/2;
zs=nz/2;

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
% coeff=[ 1.25219, -0.121591, 0.0331761, -0.0109247, 0.0035117, -0.000949881, 0.000177978];
coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordTxx=zeros(nt,nx);
seis_recordTxz=zeros(nt,nx);

% coeff2=[13/12 -5/24 1/6 -1/24];
coeff2=[1.16667, -0.416667, 0.333333, -0.0833333];

% coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];
% coeff2=[13/12 -5/24 1/6 -1/24];

Txxx_3=p;
Txxx_2=p;
Txxx_1=p;

Txzz_3=p;
Txzz_2=p;
Txzz_1=p;

Txzx_3=p;
Txzx_2=p;
Txzx_1=p;

Tzzz_3=p;
Tzzz_2=p;
Tzzz_1=p;
%%%%%%%%%%%%%%%
Vxx_3=p;
Vxx_2=p;
Vxx_1=p;

Vxz_3=p;
Vxz_2=p;
Vxz_1=p;

Vzx_3=p;
Vzx_2=p;
Vzx_1=p;

Vzz_3=p;
Vzz_2=p;
Vzz_1=p;

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

    Txzz1=(circshift(Txz,[ 1])-circshift(Txz,[ 0]));
    Txzz2=(circshift(Txz,[ 2])-circshift(Txz,[ -1]));
    Txzz3=(circshift(Txz,[ 3])-circshift(Txz,[ -2]));
    Txzz4=(circshift(Txz,[ 4])-circshift(Txz,[ -3]));
    Txzz5=(circshift(Txz,[ 5])-circshift(Txz,[ -4]));
    Txzz6=(circshift(Txz,[ 6])-circshift(Txz,[ -5]));
    Txzz7=(circshift(Txz,[ 7])-circshift(Txz,[ -6]));

    Txzz=coeff(1)*Txzz1+coeff(2)*Txzz2+coeff(3)*Txzz3+coeff(4)*Txzz4+coeff(5)*Txzz5+coeff(6)*Txzz6+coeff(7)*Txzz7;

    Tzzz1=(circshift(Tzz,[ 0])-circshift(Tzz,[ -1]));
    Tzzz2=(circshift(Tzz,[ 1])-circshift(Tzz,[ -2]));
    Tzzz3=(circshift(Tzz,[ 2])-circshift(Tzz,[ -3]));
    Tzzz4= (circshift(Tzz,[ 3])-circshift(Tzz,[ -4]));
    Tzzz5=(circshift(Tzz,[ 4])-circshift(Tzz,[ -5]));
    Tzzz6= (circshift(Tzz,[ 5])-circshift(Tzz,[ -6]));
    Tzzz7= (circshift(Tzz,[ 6])-circshift(Tzz,[ -7]));

    Tzzz=coeff(1)*Tzzz1+coeff(2)*Tzzz2+coeff(3)*Tzzz3+coeff(4)*Tzzz4+coeff(5)*Tzzz5+coeff(6)*Tzzz6+coeff(7)*Tzzz7;

    Txzx1=(circshift(Txz,[0 -1])-Txz);
    Txzx2=(circshift(Txz,[0 -2])-circshift(Txz,[0 1]));
    Txzx3=(circshift(Txz,[0 -3])-circshift(Txz,[0 2]));
    Txzx4=(circshift(Txz,[0 -4])-circshift(Txz,[0 3]));
    Txzx5=(circshift(Txz,[0 -5])-circshift(Txz,[0 4]));
    Txzx6=(circshift(Txz,[0 -6])-circshift(Txz,[0 5]));
    Txzx7=(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));

    Txzx=coeff(1)*Txzx1+coeff(2)*Txzx2+coeff(3)*Txzx3+coeff(4)*Txzx4+coeff(5)*Txzx5+coeff(6)*Txzx6+coeff(7)*Txzx7;

    Txxxabs=coeff2(1)* Txxx+coeff2(2)* Txxx_1+coeff2(3)*Txxx_2+coeff2(4)*Txxx_3;
    Tzzzabs=coeff2(1)* Tzzz+coeff2(2)* Tzzz_1+coeff2(3)*Tzzz_2+coeff2(4)*Tzzz_3;
    Txzzabs=coeff2(1)* Txzz+coeff2(2)* Txzz_1+coeff2(3)*Txzz_2+coeff2(4)*Txzz_3;
    Txzxabs=coeff2(1)* Txzx+coeff2(2)* Txzx_1+coeff2(3)*Txzx_2+coeff2(4)*Txzx_3;

    Vx=Vx+1./(rou).*dt.*(Txxxabs+Txzzabs)/h;
    Vz=Vz+1./(rou).*dt.*(Tzzzabs+Txzxabs)/h;

    Txxx_3=Txxx_2;
    Txxx_2=Txxx_1;
    Txxx_1=Txxx;


    Txzx_3=Txzx_2;
    Txzx_2=Txzx_1;
    Txzx_1=Txzx;

    Tzzz_3=Tzzz_2;
    Tzzz_2=Tzzz_1;
    Tzzz_1=Tzzz;

    Txzz_3=Txzz_2;
    Txzz_2=Txzz_1;
    Txzz_1=Txzz;
    % Vx=Vx+1./(rou).*dt.*(Txxx+Txzz)/h;
    % Vz=Vz+1./(rou).*dt.*(Tzzz+Txzx)/h;

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
save('figure8c.mat')