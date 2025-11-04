
clear;clc
close all



load('figure2_ABS_Balanced_1D_r19.mat')


% figure;plot((1:nt)*dt,analytical2,'b','LineWidth',1.5)
% hold on;plot((2:nt-1)*dt,seis_record(2:nt-1),'r','LineWidth',1.5)

figure;plot((1:nt-1)*dt,analytical2(2:end),'b','LineWidth',1.5)
hold on;plot((1:nt)*dt,seis_record,'r','LineWidth',1.5)

% load('figure1_ABS_Non_Balanced_1D_r45.mat')
% 
% 
% hold on;plot(seis_record,'b','LineWidth',1.5)


load('figure2_ABS_Non_Balanced_1D_r19.mat')



% hold on;plot((1:nt)*dt,analytical2,'r','LineWidth',1.5)
hold on;plot((1:nt)*dt,seis_record,'k--','LineWidth',1.5)

legend('Analytical','ABS\_M3','NB\_ABS\_M3')

axis([3.55 3.9 -0.5 1])
grid on 
box on
xlabel("Time in ms")
ylabel("Amplitude")



load('figure2_ABS_Balanced_1D_r38.mat')


figure;plot((1:nt)*dt,analytical2,'b','LineWidth',1.5)
hold on;plot((2:nt-1)*dt,seis_record(2:nt-1),'r','LineWidth',1.5)

% load('figure1_ABS_Non_Balanced_1D_r45.mat')
% 
% 
% hold on;plot(seis_record,'b','LineWidth',1.5)


load('figure2_ABS_Non_Balanced_1D_r38.mat')



% hold on;plot((1:nt)*dt,analytical2,'r','LineWidth',1.5)
hold on;plot((1:nt)*dt,seis_record,'k--','LineWidth',1.5)

legend('Analytical','ABS\_M3','NB\_ABS\_M3')

axis([3.55 3.9 -0.5 1])
grid on 
box on
xlabel("Time in ms")
ylabel("Amplitude")