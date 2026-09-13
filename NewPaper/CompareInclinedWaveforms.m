function CompareInclinedWaveforms(th0_list)

l=2;
m=1;

inputDB.chi1 = [0 0 0];
inputDB.geodesics = 1;
inputDB.ICs = 'spherical';
inputDB.Tmax = 1000;
inputDB.r0 = 7;
inputDB.dt = 0.25/6;
%inputDB.dt = 1;

L = length(th0_list);

for i=1:L
   th0 = th0_list(i);
   inputDB.th0 = deg2rad(th0);
   DB = DB_class(inputDB);
   hlm = DB_anal_wave(DB,2,1,'geod','rhoresum','hybrid','deltaresum','hybrid','noplot','pade');
   if i==1
       hlm_mat = repmat(hlm,1,L);
   else
       hlm_mat(:,i) = hlm;
   end
end
T = DB.t

my_linewidth = 2.0;
axes_fontsize = 12;
legend_fontsize = axes_fontsize+2;
labels_fontsize = 20;;

f = figure;
f.Position(3:4) = f.Position(3:4)*1.2;
f.Position(3) = f.Position(3).*2;
t_in = T(end)-175;
t_end = T(end)-30;
label = KerrLabel(DB);
wf_label = sprintf('h_{%d%d}',l,m);

ax1 = gca;
set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.9,'NumColumns',2)%,'Orientation','horizontal')
xlim([t_in,t_end])
ylabel(sprintf('$\\Re[%s]/\\nu$',wf_label),'Interpreter','latex','Fontsize',labels_fontsize)
hold on

for i=1:L
    th0 = th0_list(i);
    hlm_EOB = hlm_mat(:,i);
    my_linewidth = 2.15 - i.*0.15;
    plot(T,real(hlm_EOB),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\theta=%d^\\circ$',th0))
end


xs1 = xlim;
ys1 = ylim;
if abs(ys1(2)+ys1(1))>1e-3
    ys1(2) = min(abs(ys1));
    ys1(1) = -ys1(2);
end
ys1=[-m,m];

%set(ax1,'XTickLabel',[])
%text(xs1(2)-abs(xs1(2)-xs1(1))./3,ys1(2) - abs(ys1(2)-ys1(1))./10,label,'Interpreter','latex','FontSize',16)

ax2 = axes('Position',[0.6,0.20,0.25,0.40]);
set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','LineWidth',0.5,'box','on');
xlim([923,938])
hold on
for i=1:L
    th0 = th0_list(i);
    hlm_EOB = hlm_mat(:,i);
    my_linewidth = 2.15 - i.*0.15;
    plot(ax2,T,real(hlm_EOB),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\theta=%d^\\circ$',th0))
end
%for i=1:length(configs)
%    s = configs(i).s;
%    label = KerrLabel(s.dyn);
%    tLR = tLR_splined(s.dyn);
%
%    [a,t] = DB_asymmetry(s,'nodyn');
%    plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
%end
%xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',1.5*my_linewidth,'HandleVisibility','off')
%xlim(ax1,[2200,3011])
%xlim(ax2,[2200,3011])
