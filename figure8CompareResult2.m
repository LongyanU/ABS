
clear;clc;close all


load('figure8a.mat')

plotimage((1:nt-2),seis_recordVx(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordVxa=seis_recordVx;



load('figure8b.mat')


plotimage((1:nt-2),seis_recordVx(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordVxb=seis_recordVx;


load('figure8c.mat')


plotimage((1:nt-2),seis_recordVx(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordVxc=seis_recordVx;

load('figure8d.mat')
seis_recordVxd=seis_recordVx;

load('figure8e.mat')
seis_recordVxe=seis_recordVx;

figure;plot((1:nt-2),seis_recordVxa((1:nt-2),100)./max(abs(seis_recordVxa((1:nt-2),100))),'b','linewidth',2);
hold on;plot((1:nt-2),seis_recordVxb((1:nt-2),100)./max(abs(seis_recordVxb((1:nt-2),100))),'k','linewidth',2);

hold on;plot((1:nt-2),seis_recordVxc((1:nt-2),100)./max(abs(seis_recordVxc((1:nt-2),100))),'r--','linewidth',2);

hold on;plot((1:nt-2),seis_recordVxd((1:nt-2),100)./max(abs(seis_recordVxd((1:nt-2),100))),'c','linewidth',2);

hold on;plot((1:nt-2),seis_recordVxe((1:nt-2),100)./max(abs(seis_recordVxe((1:nt-2),100))),'m--','linewidth',2);


load('figure8f.mat')

% NT=nt;
% B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
% I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix
% 
% temp=I*seis_recordVx(:,100);
% hold on; plot(1:nt-2,temp(1:nt-2)./max(temp),'y','linewidth',2);


NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_recordVx(:,100);
hold on;plot((1:nt-2),temp(1:nt-2)./max(abs(temp)),'y--','linewidth',2)

legend('2-SGFD','ABSM4 Tra','NB-ABS M4','ABSM3 Tra','NB-ABS M3','TDE')
grid on

xlabel('time(ms)')
ylabel('Amp')

axis([0 998 -0.85 1.02])