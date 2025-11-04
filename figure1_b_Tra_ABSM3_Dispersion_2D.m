clear;clc;
% close all

d= [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

kh=linspace(1/1000,pi,100);

% a=[1.08333, -0.166667, 0.0833333]
a=[25/24, -2/24, 1/24]
r=0.4;
cdcfinal=zeros(5,100);
D=0;


for ii=1:100

    for aa=1:5

        xita=(aa-1)*pi/16;
        tempyy=100;
        for cdc=0.8:0.0001:1.2
            tempD1=0;
            for n=1:7
                % D=D+d(n)*sin((n-0.5)*kh(ii));
                tempD1=tempD1+d(n)* sin((n-0.5)*kh(ii)*sin(xita));
            end

            tempD2=0;
            for n=1:7
                % D=D+d(n)*sin((n-0.5)*kh(ii));
                tempD2=tempD2+d(n)* sin((n-0.5)*kh(ii)*cos(xita));
            end

            D=sqrt(tempD1^2+tempD2^2);
            temp=a(1);
            temp=temp+a(2)*cos(kh(ii)*r*cdc);
            temp=temp+a(3)*cos(2*kh(ii)*r*cdc);
            temp1=2*asin(r*D.*temp)-cdc*kh(ii)*r;
            % disp(temp1)
            if abs(temp1)<tempyy
                tempyy=abs(temp1);
                cdcfinal(aa,ii)=cdc;
            end
        end
    end
end


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
% ylabel('c_fd/c')
ylabel('$\frac{c_{\mathrm{fd}}}{c}$', 'Interpreter', 'latex');
xlim([0,pi])
% ylim([0.89,1.02])
ylim([0.94,1.08])