clear;
close all;
clc

% clip level 0.5


% load('figure5e.mat')
% plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')
% 
% 
% NT=nt;
% B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
% I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix
% 
% temp=I*seis_record(:,75);



load('figure5a.mat')
seis_recorda=seis_record;
plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')


load('figure5b.mat')
seis_recordb=seis_record;
% plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')



load('figure5c.mat')
% plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')

seis_recordc=seis_record;


load('figure5d.mat')
% plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')

seis_recordd=seis_record;


load('figure5e.mat')
% plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')

seis_recorde=seis_record;

figure;plot((1:nt-2)*1,seis_recorda((1:nt-2),75),'b','linewidth',2)

hold on;plot((1:nt-2)*1,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*1,seis_recordc((1:nt-2),75),'k--','linewidth',2)

hold on;plot((1:nt-2)*1,seis_recordd((1:nt-2),75),'c','linewidth',2)
hold on;plot((1:nt-2)*1,seis_recorde((1:nt-2),75),'m--','linewidth',2)



load('figure5f.mat')
% plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')


NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);

hold on; plot((1:nt-2)*1,temp((1:nt-2)),'y:','linewidth',2);


xlabel('time(ms)')
ylabel('Amp')
legend('2-SGFD','ABSM4 Tra','NB-ABS M4','ABSM3 Tra','NB-ABS M3','TDE')
grid on
axis([ 0 5300  -1.15*10^-5 1.95*10^-5])
box on

% axis([ 4900 5050  -0.012 0.02])
% 
% axis([ 4300 4500  -0.05 0.05])