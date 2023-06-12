% population model with spatial heterogenous chi
clear;
%% parameters for ligand-receptor bindings 
KI1 = 18.2;		%Dissociation constant between inactive Tar receptor and MeAsp, unit:um
KA1 = 3000;
KI2 = 70;		%Dissociation constant between inactive Tar receptor and MeAsp, unit:um
KA2 = 1000;
KI3 = 1500;		%Dissociation constant between inactive Tar receptor and MeAsp, unit:um
KA3 = 3500;

N=6;            %receptor number in an effective cluster
coff_ita0=15;    
%% attractant gradients
dx=0.88;   x=(1:850)*dx;                      % channel length 
attr=2000; L = (1-x/max(x))*attr/3 + attr/3;  % MeAsp gradient
attr2=2.^(0:1:18);                            % loading concentration of serine

%% parameters for lognormal distribution of Tar/Tsr
mu = log(.25); 
sigma = 1.5; 

%% cell distribution of steady state;
Rho=zeros(length(attr2),length(x));
veff=zeros(length(attr2),length(x));

for ii = 1:length(attr2)
L2 = x/max(x)*attr2(ii)/3+attr2(ii)/3;        % serine gradient
N2=1:100;
ratio = 0.5;
fL=(N* ratio./(1+ratio))*log((1+L/KI1) ./(1+L/KA1));
fL2=(N./(1+ratio))*log((1+L2/KI2)./(1+L2/KA2));
fL3=(N* ratio./(1+ratio))*log((1+L2/KI3)./(1+L2/KA3));
fL_Total=fL+fL2+fL3;

% % consider diversity of tumble bias
% bias=0:0.1:0.5;

% fis bias=0.2
bias = 0.2;
rho_temp=zeros(length(bias),length(x));
rho=zeros(1,length(x));
cmap2=winter(length(bias));

for iii=1:length(bias)
    coff_ita=coff_ita0*(1-bias(iii));
    ro=exp(coff_ita.*fL_Total);
    rho_temp(iii,:)=ro/sum(ro);
    P3=normpdf(bias(iii),0.2,0.05);               %% lognormal 
    rho=rho+P3*rho_temp(iii,:);
    rho=rho./sum(rho);
% switch j
%     case 1
%        plot(x,rho_temp(iii,:),'color',[1,(iii)/(1+length(bias)),(iii)/(1+length(bias))],'linewidth',2); hold on; 
%     case 2
%        plot(x,rho_temp(iii,:),'color',[(iii)/(1+length(bias)),1,(iii)/(1+length(bias))],'linewidth',2); hold on;
%     case 3
%        plot(x,rho_temp(iii,:),'color',[(iii)/(1+length(bias)),(iii)/(1+length(bias)),1],'linewidth',2); hold on;
% end
end
% Rho(ii,:) = Rho(ii,:)./sum(Rho(ii,:));
Rho(ii,:) = rho;
end

%% plot heat map
attr2=[attr2,attr2(end)];
Rho = [Rho;Rho(end,:)];

figure('color','w');
surf(x,attr2,Rho); shading flat; view(2);
xlabel('x (\mum)'); xlim([0,748]); xticks([0,200,400,600]);
ylabel('Serine Conc. (\muM)'); ylim([4,64000]); yticks([10,100,1000,10000]);
set(gcf,'unit','centimeters','position',[10 10 8 7]);
set(gca,'yscale','log','FontName','Arial','FontSize',9);
set(gca,'Position',[0.14,0.15,0.7,0.75]);
colormap jet;  h=colorbar;  caxis([0,0.01]);
ax=gca; axpos=ax.Position; h.Position(3)=.7*h.Position(3); h.Position(1)=0.86;
ax.Position=axpos; h.Ticks=0:0.002:0.01; h.Ruler.Exponent=-3;

% figure('color','w');  plot(x,Rho([4,7,11,14],:)','linewidth',2); ylim([0 0.01]);  xlim([0 740]); xlabel('x (\mum)'); ylabel('Pdf');  


%% calculate w1,w2
x=x';
Me=mean(Rho(1:end-1,:),2); 
for i=1:size(Rho,1)-1
    ind=find(Rho(i,:)>Me(i)*0.9);
    Ex1(i)=sum(x(ind).*Rho(i,ind)')/max(x)./sum(Rho(i,ind));
    Ex2(i)=sqrt(sum((x(ind)-max(x)/2).^2.*(Rho(i,ind)'-Me(i)*0.9))./max(x).^2); 
end

figure('color','w'); subplot(2,1,1); 
semilogx(attr2(3:end-1),Ex1(3:end),'o-','color','k','markersize',5,'markeredgecolor','k'); 
line(2.^(2:18),[0.3*ones(1,17);0.7*ones(1,17)],'color',[0.5,.5,.5],'linestyle','-'); 
ylim([0,1]); xlim([5,50000]); ylabel('w_1');
hold on; subplot(2,1,2);
semilogx(attr2(3:end-1),Ex2(3:end),'o-','color','k','markersize',5,'markeredgecolor','k'); 
line(2.^(2:18),0.15*ones(1,17),'color',[0.5,.5,.5],'linestyle','-'); 
ylim([0,0.4]);  ylabel('w_2'); 
xlim([5,50000]); xlabel('Serine conc. (\muM)');

set(gcf,'unit','centimeters','position',[10 5 7.5 8]);
