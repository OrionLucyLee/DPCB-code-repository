%% load experimental data
clc; clear;

x_1 = xlsread('RatioDistribution.csv','A2:A41');
Dis_1 = xlsread('RatioDistribution.csv','B2:B41');
lgx_1 = log(x_1);
Rholog_1 = x_1.*Dis_1*10;

x_2 = xlsread('RatioDistribution.csv','A43:A82');
Dis_2 = xlsread('RatioDistribution.csv','B43:B82');
lgx_2 = log(x_2);
Rholog_2 = x_2.*Dis_2*10;

x_3 = xlsread('RatioDistribution.csv','A84:A123');
Dis_3 = xlsread('RatioDistribution.csv','B84:B123');
lgx_3 = log(x_3);
Rholog_3 = x_3.*Dis_3*10;

% figure; plot(x_1,Dis_1,'or',x_2,Dis_2,'ob',x_3,Dis_3,'og'); figure;
% plot(lgx_1,Rholog_1,'or',lgx_2,Rholog_2,'ob',lgx_3,Rholog_3,'og');

%% gaussian fitting
x=-5:.1:3;
Mytype = fittype('A*exp(-(x-u)^2/(2*d^2))');
[cf1,gof1] = fit(lgx_1,Rholog_1,Mytype);
y1 = cf1.A.*exp(-(x-cf1.u).^2/(2*cf1.d.^2));
[cf2,gof2] = fit(lgx_2,Rholog_2,Mytype);
y2 = cf2.A.*exp(-(x-cf2.u).^2/(2*cf2.d.^2));
[cf3,gof3] = fit(lgx_3,Rholog_3,Mytype);
y3 = cf3.A.*exp(-(x-cf3.u).^2/(2*cf3.d.^2));

% three separated gaussian distributions
figure('color','w'); set(gcf,'unit','centimeters','position',[10 5 8 7]);
scatter(lgx_1,Rholog_1,15,'LineWidth',.5,'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none', 'MarkerFaceColor','m'); hold on;
plot(x,y1,'r','linewidth',2); 
scatter(lgx_2,Rholog_2,15,'LineWidth',.5, 'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none','MarkerFaceColor','g'); hold on;
plot(x,y2,'g','linewidth',2);
scatter(lgx_3,Rholog_3,15,'LineWidth',.5, 'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none','MarkerFaceColor','b'); hold on;
plot(x,y3,'b','linewidth',2);
xlim([-4.5,2.5]); set(gca,'Fontsize',10); box on;
xlabel('log(r)'); ylabel('Probability');    

% scaled distribution
figure('color','w');  set(gcf,'unit','centimeters','position',[10 5 4,3.5]);
plot((x-cf1.u)/cf1.d, y1/cf1.A, 'linewidth',1);
hold on; scatter((lgx_1-cf1.u)/cf1.d, Rholog_1/cf1.A,5,'LineWidth',.5,'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none', 'MarkerFaceColor','r'); hold on;
scatter((lgx_2-cf2.u)/cf2.d, Rholog_2/cf2.A,5,'LineWidth',.5,'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none', 'MarkerFaceColor','g');
scatter((lgx_3-cf3.u)/cf3.d, Rholog_3/cf3.A,5,'LineWidth',.5,'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none', 'MarkerFaceColor','b');
xlim([-4,4]); set(gca,'Fontsize',8); box on;
xlabel('(log(r)-\mu)/\sigma'); ylabel('Scaled Probability');    

%% fitting of mu and sigma
od=[0.05,0.3,0.51];
A = [cf1.A,cf2.A,cf3.A];
A_lower = [ 0.4533, 0.6302, 0.8142];
A_upper = [ 0.4949, 0.6692, 0.8472];

mu = [cf1.u,cf2.u,cf3.u];
mu_lower = [ -0.7443, 0.02187, 0.4152];
mu_upper = [ -0.6515, 0.06535, 0.4372];

sigma = abs([cf1.d,cf2.d,cf3.d]);
sigma_lower = [0.7147,  0.5373, 0.4327];
sigma_upper = [0.79,    0.577,  0.4544];

t_A = polyfit(od,A,1);
t_mu = polyfit(od,mu,1);
t_sigma = polyfit(od,sigma,1);

figure('color','w'); set(gcf,'unit','centimeters','position',[10 5 4,3.5]);
errorbar(od,A, A_lower-A, A_upper-A,'or','linewidth',1.25,'markersize',4);
hold on; plot(od,polyval(t_A,od),'k','linewidth',1.5);
xlim([0,0.6]); ylim([0.03,0.1]); set(gca,'Fontsize',10); box on;
xlabel('OD_6_0_0'); ylabel('A');         

figure('color','w'); set(gcf,'unit','centimeters','position',[10 5 4.5,8]);
subplot(2,1,1,'Position',[0.1,0,.8,.4]);
errorbar(od,mu, mu_lower-mu, mu_upper-mu,'or','linewidth',1.25,'markersize',4);
hold on; plot(od,polyval(t_mu,od),'k','linewidth',1.5);
xlim([0,0.6]); ylim([-1.1,1]); ylabel('\mu');  set(gca,'Fontsize',10,'xtick',[]); box on;
        
subplot(2,1,2,'Position',[0.1,0.5,.8,.4]);
errorbar(od,sigma, sigma_lower-sigma, sigma_upper-sigma,'or','linewidth',1.25,'markersize',4);
hold on; plot(od,polyval(t_sigma,od),'k','linewidth',1.5);
xlim([0,0.6]); ylim([0.35, .85]); set(gca,'Fontsize',10); box on;
xlabel('OD_6_0_0'); ylabel('\sigma');   

%% hypothetic distribution for OD=0.05,0.2,0.3
save('FittingPara','od','t_A','t_mu','t_sigma');
od_test=[0.05,0.2,0.35];
AA=polyval(t_A,od_test); uu=polyval(t_mu,od_test);dd=polyval(t_sigma,od_test); 
y1=AA(1).*exp(-(x-uu(1)).^2/(2*dd(1).^2)); 
y2=AA(2).*exp(-(x-uu(2)).^2/(2*dd(2).^2)); 
y3=AA(3).*exp(-(x-uu(3)).^2/(2*dd(3).^2)); 
dis_mix = (y1+y3)/2;
figure('color','w'); set(gcf,'unit','centimeters','position',[10 5 8 7]);
plot(x,y1,x,y2,x,y3, x,dis_mix);

% average and variance
m13=sum(x.*dis_mix)/sum(dis_mix); 
st13=sqrt(sum(dis_mix.*(x-m13).^2)/(sum(dis_mix)));

figure('color','w'); set(gcf,'unit','centimeters','position',[10 5 8 7]);
plot(od_test,[uu(1),uu(2),uu(3)],'o-'); hold on; plot(0.2,m13,'o'); 
set(gca,'Fontsize',10); xlabel('OD_6_0_0'); ylabel('Average log(r)');
ylim([-1.1,0.3]);

figure('color','w'); set(gcf,'unit','centimeters','position',[10 5 8 7]);
plot(od_test,[dd(1),dd(2),dd(3)],'v-'); hold on; plot(0.2,st13,'v'); 
set(gca,'Fontsize',10); xlabel('OD_6_0_0'); ylabel('Variance of log(r)');


