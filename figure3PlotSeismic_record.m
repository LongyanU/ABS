
clear;clc;close all

load('figure3_tra.mat')
plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')


load('figure3b.mat')
plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')


load('figure3c.mat')
plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

load('figure3d.mat')
plotimage((1:nt-2)*1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')