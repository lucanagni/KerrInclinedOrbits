function DB_multipolar(s)
    t = s.ell(2).emm(3).t;
    
    figure
    tiledlayout(1,2)
    
    nexttile
    hold on
    plot(t,abs(s.ell(2).emm(3).hlm),'black','LineStyle','-','DisplayName','(2,2)')
    plot(t,abs(s.ell(2).emm(2).hlm),'black','LineStyle','--','DisplayName','(2,1)')
    plot(t,abs(s.ell(3).emm(4).hlm),'blue','LineStyle','-','DisplayName','(3,3)')
    plot(t,abs(s.ell(3).emm(3).hlm),'blue','LineStyle','--','DisplayName','(3,2)')
    plot(t,abs(s.ell(4).emm(5).hlm),'red','LineStyle','-','DisplayName','(4,4)')
    legend
    
    xlim([3850,4100])
    xlabel('t')
    ylabel('|h_{lm}|')


    nexttile
    hold on
    plot(t,freq(s.ell(2).emm(3).hlm,t),'black','LineStyle','-','DisplayName','(2,2)')
    plot(t,freq(s.ell(2).emm(2).hlm,t),'black','LineStyle','--','DisplayName','(2,1)')
    plot(t,freq(s.ell(3).emm(4).hlm,t),'blue','LineStyle','-','DisplayName','(3,3)')
    plot(t,freq(s.ell(3).emm(3).hlm,t),'blue','LineStyle','--','DisplayName','(3,2)')
    plot(t,freq(s.ell(4).emm(5).hlm,t),'red','LineStyle','-','DisplayName','(4,4)')
    legend

    xlim([3850,4100])
    xlabel('t')
    ylabel('\omega_{lm}')

return