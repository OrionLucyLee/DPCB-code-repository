% code for simulating the transitions curves 

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
mu = -0.7; 
sigma = .75;

%% plot the log-normal distribution
xx = exp(-5:.1:5);
pdf = lognpdf(xx, mu,sigma);
figure('color','w'); 
set(gcf,'unit','centimeters','position',[10 5 4,3]);
plot(log(xx),xx.*pdf);
xlim([-4,2.5]); set(gca,'Fontsize',9); box on;
xlabel('ln(Tar/Tsr)'); ylabel('Probability');    


%% vary f_c and plot the two transitin curves
f_c = [0.01,0.05,0.1, 0.2, 0.3];

r1=logninv(f_c,mu,sigma);
figure;set(gcf,'unit','centimeters','position',[10 5 8 7]);
for i=1:length(r1)
   yy=r1(i)*cm/2.*a1./(a2+r1(i).*a3);
   xx=z-yy; 
   xx(find(xx<0))=0.1; 
   loglog(xx,yy,'color',[1-i/6,1-i/6,1],'linewidth',2);
   hold on; 
  legend_str{i}= num2str(f_c(i));
end
legend(legend_str);
xlabel('c_m_i_n (\muM)','FontName','Arial','FontSize',10); xlim([4,30000]); xticks([10,100,1000,10000]); 
ylabel('c_m_a_x-c_m_i_n (\muM)','FontName','Arial','FontSize',10); ylim([4,30000]); box on; yticks([10,100,1000,10000]); 


r2=logninv(1-f_c,mu,sigma);
figure;set(gcf,'unit','centimeters','position',[10 5 8 7]);
for i=1:length(r2)
   y=r2(i)*cm/2.*A1./(A2+r2(i).*A3);
   loglog(x,y,'color',[1,1-i/6,1-i/6],'linewidth',2);
   hold on; 
   legend_str2{i}= num2str(f_c(i));
end
legend(legend_str2);
xlabel('c_m_i_n (\muM)','FontName','Arial','FontSize',10); xlim([4,30000]); xticks([1,10,100,1000,10000]); 
ylabel('c_m_a_x-c_m_i_n (\muM)','FontName','Arial','FontSize',10); ylim([4,30000]); box on; yticks([1,10,100,1000,10000]); 



