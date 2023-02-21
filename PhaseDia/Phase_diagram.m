% code for obtaining the phase diagram 

clc; clear;
%% parameters for ligand-receptor bindings 
KI1 = 18.2;		%Dissociation constant of MeAsp binding to inactive Tar, unit:uM
KA1 = 3000;     %Dissociation constant of MeAsp binding to active Tar, unit:uM
KI2 = 70;		%Dissociation constant of serine binding to inactive Tsr, unit:uM
KA2 = 1000;     %Dissociation constant of serine binding to active Tsr, unit:uM
KI3 = 1500;		%Dissociation constant of MeAsp binding to inactive Tsr, unit:uM
KA3 = 3500;     %Dissociation constant of MeAsp binding to active Tsr, unit:uM

coff_ita=10; 
cm=2000;       % loading conc. of MeAsp
gm=cm/3;       % drop of MeAsp gradient


%% parameters for log-normal distribution of Tar/Tsr
f_c= 0.1;   % f_c
mu = -0.843; 
sigma = .97; 


%% plot the log-normal distribution
xx = exp(-5:.1:5);
pdf = lognpdf(xx, mu,sigma);
figure('color','w'); 
set(gcf,'unit','centimeters','position',[10 5 4,3]);
plot(log(xx),xx.*pdf);
xlim([-4,2.5]); set(gca,'Fontsize',9); box on;
xlabel('ln(Tar/Tsr)'); ylabel('Probability');    


%%  phase diagram
r0=logninv(1-f_c,mu,sigma);
rL=logninv(f_c,mu,sigma);

% r_min = r_1
x=2.^(0:15);   % c_min
A1=(KA1-KI1)./((cm*2/3+KI1).*(cm*2/3+KA1));
A2=(KA2-KI2)./((x+KI2).*(x+KA2));
A3=(KA3-KI3)./((x+KI3).*(x+KA3));
y=r0*gm.*A1./(A2+r0.*A3);

% r_max = r_2
z=2.^(2:0.01:15);   % c_max
a1=(KA1-KI1)./((cm/3+KI1).*(cm/3+KA1));
a2=(KA2-KI2)./((z+KI2).*(z+KA2));
a3=(KA3-KI3)./((z+KI3).*(z+KA3));
yy=rL.*gm.*a1./(a2+rL.*a3);
xx=z-yy; 
xx(find(xx<0))=0.1; 

figure('color','w'); set(gcf,'unit','centimeters','position',[10 10 8 7]);
loglog(x,y,'color',[1,0.2,0.2],'linewidth',2); hold on;
loglog(xx,yy,'color',[0.2,0.2,1],'linewidth',2);
xlabel('c_m_i_n (\muM)','FontName','Arial','FontSize',12); xlim([4,30000]); xticks([10,100,1000,10000]); 
ylabel('c_m_a_x-c_m_i_n (\muM)','FontName','Arial','FontSize',12); ylim([4,30000]); box on; yticks([10,100,1000,10000]); 
ax=gca; load('MyColormaps.mat'); colormap(ax,mycmap); set(gca,'xscale','log','yscale','log','FontName','Arial','FontSize',10);


%% chemotaxis states from experiments with cells collected at OD=0.1
S1_x = [4,8,16,12000,16000]/3;  
S2_x = [32,64,128,256]/3;
S3_x = [512,1000,1500,2000,3000,4000]/3;  
S4_x = 8000/3;
scatter(S1_x,S1_x,20,'MarkerEdgeColor','k', 'MarkerFaceColor','r', 'LineWidth',1); hold on;
scatter(S2_x,S2_x,20,'MarkerEdgeColor','k', 'MarkerFaceColor','y', 'LineWidth',1); hold on;
scatter(S3_x,S3_x,20,'MarkerEdgeColor','k', 'MarkerFaceColor','b', 'LineWidth',1); hold on;
scatter(S4_x,S4_x,20,'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'LineWidth',1); hold off;
