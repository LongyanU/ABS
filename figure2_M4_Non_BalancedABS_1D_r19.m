% 历时 38.682071 秒。
% 历时 38.311698 秒。
% 历时 37.785946 秒。
% 历时 33.637174 秒。
clear
clc %%%%%%%
close all
% Elapsed time is 9.253353 seconds.
nt=2400;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling

nx=850;


v=2000;
c=v;
% v(1:nz/2,:)=2800;


p=zeros([1,nx]); pboundarynew=p;pdan=p;
dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;

dt=0.0019; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0 = 9;
t0 = 4/(2*pi*f0);       % 信号峰值时间 (s)
t = (0:nt-1)*dt;        % 时间轴 t ∈ [0, (nt-1)*dt]

src = (1 - 2*(pi*f0*(t - t0 )).^2) .* ...
    exp(-(pi*f0*(t - t0 )).^2);


seis_record=zeros(nt,1);
p=zeros([nx,1]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
xs=nx/2;

h=dx;
r=v*dt/h;

d1pxz_2=p;
d1pxz_3=p;
d1pxz_4=p;





% coeff2=[1.15972, -0.395833, 0.3125, -0.0763889] 
coeff2=[1.16667, -0.416667, 0.333333, -0.0833333]

coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];
% figure;
for loop=1:1
    tic

    p=zeros([nx,1]); Vx=p; Vz=p;

    for it=1:nt,

        d1px11=Vx-circshift(Vx,[1]);
        d1px12=(circshift(Vx,[-1])-circshift(Vx,[2]));
        d1px13=(circshift(Vx,[-2])-circshift(Vx,[3]));
        d1px14=(circshift(Vx,[-3])-circshift(Vx,[4]));
        d1px15=(circshift(Vx,[-4])-circshift(Vx,[5]));
        d1px16=(circshift(Vx,[-5])-circshift(Vx,[6]));
        d1px17=(circshift(Vx,[-6])-circshift(Vx,[7]));

        d1px=coeff(1)*d1px11+coeff(2)*d1px12+coeff(3)*d1px13+coeff(4)*d1px14+coeff(5)*d1px15+coeff(6)*d1px16...
            +coeff(7)*d1px17;


        d1pxz=d1px;
        p=p-dt*v.^2.*(coeff2(1)*d1pxz+coeff2(2)*d1pxz_2+coeff2(3)*d1pxz_3+coeff2(4)*d1pxz_4)/h;

        d1pxz_4=d1pxz_3;
        d1pxz_3=d1pxz_2;
        d1pxz_2=d1pxz;

        p(xs)= src(it);


        d1px1=(circshift(p,[-1])-circshift(p,[0]));
        d1px2=(circshift(p,[-2])-circshift(p,[1]));
        d1px3=(circshift(p,[-3])-circshift(p,[2]));
        d1px4=(circshift(p,[-4])-circshift(p,[3]));
        d1px5=(circshift(p,[-5])-circshift(p,[4]));
        d1px6=(circshift(p,[-6])-circshift(p,[5]));
        d1px7=(circshift(p,[-7])-circshift(p,[6]));


        d1px=coeff(1)*d1px1+coeff(2)*d1px2+coeff(3)*d1px3+coeff(4)*d1px4+coeff(5)*d1px5+coeff(6)*d1px6...
            +coeff(7)*d1px7;


        Vx=Vx-dt*(d1px)/h;


        if rem(it,isnap)== 0,
            plot(p)
            title(it)
            drawnow
        end

        seis_record(it)=p(60);


    end

end
toc




x = 3650*2;               % 震源-接收器距离 (m)

%% 计算解析解
delta_t = x / c;        % 时间延迟
analytical = (1/(2*c)) * (1 - 2*(pi*f0*(t - t0 - delta_t)).^2) .* ...
    exp(-(pi*f0*(t - t0 - delta_t)).^2);

%% 归一化（匹配数值导数振幅）
analytical2 = analytical / max(abs(analytical)) * max(abs(src));
figure;plot(analytical2,'r')
hold on;plot(seis_record,'k')
legend('Analytical2','ABS')

% coeff2=[1.13889, -0.333333, 0.25, -0.0555556];
% save('figure1_ABS_Non_Balanced_1D.mat')


% [1.125, -0.291667, 0.208333, -0.0416667]
% save('figure1_ABS_Non_Balanced_1D_1.mat')


% [1.15972, -0.395833, 0.3125, -0.0763889] 
% save('figure1_ABS_Non_Balanced_1D_2.mat')


% [0.840278, 0.145833, 0.0208333, -0.00694444]
save('figure2_M4ABS_Non_Balanced_1D_r19.mat')