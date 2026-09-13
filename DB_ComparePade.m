function DB_ComparePade(a)

%==========================================================================
% Compare approximants for spinless test paticle 
% a is the spin of the BH
%==========================================================================

x = linspace(0,0.45,1000);

figure
grid on
hold on
for i = 0:7
    m = i;
    n = 7-i;
    if m==7
        name = sprintf('Taylor');
    else
        name = sprintf('(%d,%d)',m,n);
    end
    plot(x,DB_PlotFlux(m,n,x,a),'LineWidth',3, 'DisplayName',name)
end

hold off

lgd = legend('Location','northwest');
lgd.FontSize = 20;
xlabel('$v$','Interpreter','latex');
ylabel('$\hat{f}$','Interpreter','latex');

return

function y = DB_PlotFlux(m,n,v_omg,a)

q=0;

nu         = q/(1+q).^2;
eulergamma = 0.57721566490153286061;
X1         = q./(1+q);
X2         = 1./(1+q);

chi1 = [0;0;0];
chi2 = [0;0;a];

l = [0;0;1];

f2    = -1247./336 -35./12.*nu;
f3    = 4.*pi;
f4    = -44711./9072 + 9271./504.*nu + 65./18.*nu.^2;
f5    = -(8191./672 + 583./24.*nu).*pi;
f6    = 6643739519./69854400 + 16./3.*pi.^2 - 1702./105.*eulergamma...
    -(134543./7776 - 41./48.*pi.^2).*nu - 94403./3024.*nu.^2 - 775./324.*nu.^3;
fl6   = -1712./105;
f7    = -(16285./504 - 214745./1728.*nu - 193385./3024.*nu.^2).*pi;

f3so  = -0.25.*(11.*X1 + 5.*X2).*X1.*dot(l,chi1)...
    -0.25.*(11.*X2 + 5.*X1).*X2.*dot(l,chi2);
f4ss  = nu./48.*(289.*dot(l,chi1).*dot(l,chi2) - 103.*dot(chi1,chi2));

y = DB_pade(m,n,f2,f3,f4,f5,f6,fl6,f7,f3so,f4ss,v_omg);

return