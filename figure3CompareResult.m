clear;clc;close all

load('figure3_tra.mat')
figure;plot((1:nt-2),seis_record((1:nt-2),75),'c','linewidth',2);

load('figure3b.mat')
hold on;plot((1:nt-2),seis_record((1:nt-2),75),'k','linewidth',2);


load('figure3c.mat')
hold on;plot((1:nt-2),seis_record((1:nt-2),75),'r-.','linewidth',2);

load('figure3d.mat')
hold on;plot((1:nt-2),seis_record((1:nt-2),75),'m','linewidth',2);

load('figure3e.mat')
hold on;plot((1:nt-2),seis_record((1:nt-2),75),'y--','linewidth',2);

load('figure3f.mat')
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);
hold on;plot((1:nt-2),temp((1:nt-2)),'b:','linewidth',2);

legend('2-SGFD','Tra ABS4','non-balanced in time ABS4','Tra ABS-3','NB ABS-3','TDE')
grid on
xlabel('time(ms)')
ylabel('Amp')
box on
xlim([0 3000])
ylim([-1.2*10^-5 2*10^-5])