clear;clc;
close all

d= [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

kh=linspace(1/1000,pi,100);

% a=[1.08333, -0.166667, 0.0833333]
r=0.34;
cdcfinal=zeros(5,100);
D=0;


for ii=1:100

    for aa=1:5

        xita=(aa-1)*pi/16;

        tempD1=0;
        for n=1:7
            tempD1=tempD1+d(n)* sin((n-0.5)*kh*sin(xita));
        end

        tempD2=0;
        for n=1:7
            tempD2=tempD2+d(n)* sin((n-0.5)*kh*cos(xita));
        end

        D=sqrt(tempD1.^2+tempD2.^2);

       cdcfinal(aa,:) =2/r*asin(r*D)./kh;

    end
end

  % delta=2./(r*k*h);
  %   delta=delta.*asin(r*sqrt(a+b));
  %   a1=delta;

for ii=1:5
    if ii==1
        figure;plot(kh,cdcfinal(ii,:),'r','LineWidth',1)
    elseif ii==2
        hold on;plot(kh,cdcfinal(ii,:),'b','LineWidth',1)
    elseif ii==3
        hold on;plot(kh,cdcfinal(ii,:),'c','LineWidth',1)
    elseif ii==4
        hold on;plot(kh,cdcfinal(ii,:),'k','LineWidth',1)
    else
        hold on;plot(kh,cdcfinal(ii,:),'m','LineWidth',1)
    end
end

grid on
legend('\theta=0','\theta=\pi/16','\theta=2\pi/16','\theta=3\pi/16','\theta=4\pi/16')

xlabel('kh')
ylabel('$\frac{c_{\mathrm{fd}}}{c}$', 'Interpreter', 'latex');
axis([0 100 -1.5*10^-5  5*10^-5])
xlim([0,pi])
% ylim([0.89,1.02])
ylim([0.9,1.05])