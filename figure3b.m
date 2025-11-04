% 历时 63.651948 秒。
% 历时 66.984950 秒。
clear
clc %%%%%%%
close all

nt=1500*2;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling

nx=450;
nz=450;

v=ones(nz,nx)*3400;
v(1:nz/2,:)=2800;


p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=21.5*pi;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian



seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=60;
xs=nz/2;

h=dx;
r=v*dt/h;
p=zeros([nz nx]); Vx=p; Vz=p;

d1pxz_2=p;
d1pxz_3=p;
d1pxz_4=p;


d1px_2=p;
d1pz_2=p;

d1px_3=p;
d1pz_3=p;

d1px_4=p;
d1pz_4=p;

coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348]
coeff2=[13/12 -5/24 1/6 -1/24];
tic
for it=1:nt-2,

    d1px11=Vx-circshift(Vx,[0 1]);
    d1px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
    d1px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
    d1px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
    d1px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
    d1px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
    d1px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));


    d1pz11=Vz-circshift(Vz,[1 0]);
    d1pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
    d1pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
    d1pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
    d1pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
    d1pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
    d1pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));



    d1px=coeff(1)*d1px11+coeff(2)*d1px12+coeff(3)*d1px13+coeff(4)*d1px14+coeff(5)*d1px15+coeff(6)*d1px16...
        +coeff(7)*d1px17;
    d1pz=coeff(1)*d1pz11+coeff(2)*d1pz12+coeff(3)*d1pz13+coeff(4)*d1pz14+coeff(5)*d1pz15+coeff(6)*d1pz16...
        +coeff(7)*d1pz17;


 d1pxz=d1pz+d1px;
    p=p-dt*v.^2.*(coeff2(1)*d1pxz+coeff2(2)*d1pxz_2+coeff2(3)*d1pxz_3+coeff2(4)*d1pxz_4)/h;


    d1pxz_4=d1pxz_3;
    d1pxz_3=d1pxz_2;
    d1pxz_2=d1pxz;
   
    p(zs,xs)= p(zs,xs)+src(it)*dt^2;
    % [p,p]=spongeABC(p,p,nx,nz,50,50,0.007);

    d1px1=(circshift(p,[0 -1])-circshift(p,[0 0]));
    d1px2=(circshift(p,[0 -2])-circshift(p,[0 1]));
    d1px3=(circshift(p,[0 -3])-circshift(p,[0 2]));
    d1px4=(circshift(p,[0 -4])-circshift(p,[0 3]));
    d1px5=(circshift(p,[0 -5])-circshift(p,[0 4]));
    d1px6=(circshift(p,[0 -6])-circshift(p,[0 5]));
    d1px7=(circshift(p,[0 -7])-circshift(p,[0 6]));

    d1pz1=(circshift(p,[-1])-circshift(p,[0]));
    d1pz2=(circshift(p,[-2])-circshift(p,[1]));
    d1pz3=(circshift(p,[-3])-circshift(p,[2]));
    d1pz4=(circshift(p,[-4])-circshift(p,[3]));
    d1pz5=(circshift(p,[-5])-circshift(p,[4]));
    d1pz6=(circshift(p,[-6])-circshift(p,[5]));
    d1pz7=(circshift(p,[-7])-circshift(p,[6]));


    d1px=coeff(1)*d1px1+coeff(2)*d1px2+coeff(3)*d1px3+coeff(4)*d1px4+coeff(5)*d1px5+coeff(6)*d1px6...
        +coeff(7)*d1px7;
    d1pz=coeff(1)*d1pz1+coeff(2)*d1pz2+coeff(3)*d1pz3+coeff(4)*d1pz4+coeff(5)*d1pz5+coeff(6)*d1pz6...
        +coeff(7)*d1pz7;

    Vx=Vx-dt*(coeff2(1)*d1px+coeff2(2)*d1px_2+coeff2(3)*d1px_3+coeff2(4)*d1px_4)/h;
    Vz=Vz-dt*(coeff2(1)*d1pz+coeff2(2)*d1pz_2+coeff2(3)*d1pz_3+coeff2(4)*d1pz_4)/h;

    d1px_4=d1px_3;
    d1px_3=d1px_2;
    d1px_2=d1px;

    d1pz_4=d1pz_3;
    d1pz_3=d1pz_2;
    d1pz_2=d1pz;

    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,50,50,0.007);


    if rem(it,isnap)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end


    seis_record(it,:)=p(zs,:);


end
toc
save('figure3b.mat')
% figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
% xlabel('time(ms)')
% ylabel('Amp')
% legend('receiver A','receiver B')
% grid on
