function ISCO_position(varargin)
    % =====================================================================================================================================
    % Reproduce Fig. 1 of https://arxiv.org/pdf/1901.05901
    % =====================================================================================================================================
    type = 'LSSO';
    if ~isempty(varargin)
        type = varargin{1};
    end

    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;


    A = [0,0.2,0.5,0.9];
    %A = 0:0.01:0.99;
    I = linspace(0,180,180);
    r_LSO = zeros(length(I),length(A));

    for j=1:length(A)
        for i=1:length(I)
            fprintf('a = %.1f, iota = %3.f\n',A(j),I(i))
            if strcmp(type,'LSSO')
                r_LSO(i,j) = DB_LSSO(A(j),I(i));
                Colors = autumn(length(A));
            elseif strcmp(type,'LR')
                r_LSO(i,j) = DB_LR(A(j),I(i));
                Colors = winter(length(A));
            else
                error('varargin must be either LR or LSSO')
            end
        end
    end

    f = figure;
    f.Position(4) = f.Position(3).*1.02;
    set(groot, 'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth', 1);%,'TickLength', [0.015 0.025] %font used to be Times New Roman
    %ylabel('$\cos(\theta_\mathtt{inc})$','Interpreter','latex','FontSize',labels_fontsize)
    xlabel('$\cos\iota$','Interpreter','latex','FontSize',labels_fontsize)
    xlim([-1.01,1.01])
    xticks(-1:0.5:1);
    if strcmp(type,'LSSO')    
        ylim([1.9,9.1])
        %xlabel('\textrm{r}$_\mathtt{{LSSO}}$','Interpreter','latex','FontSize',labels_fontsize)
        ylabel('r$_\mathrm{{LSSO}}$','Interpreter','latex','FontSize',labels_fontsize)
    else
        ylim([1.4,4.1])
        %xlabel('\textrm{r}$_\mathtt{{LR}}$','Interpreter','latex','FontSize',labels_fontsize)
        ylabel('r$_\mathrm{{LR}}$','Interpreter','latex','FontSize',labels_fontsize)
    end

    %Colors = ['k','g','c','r'];
    %Colors = [[0 0 0];MyColors('b1');MyColors('r1');MyColors('g1')];
    LineStyles = ["-" "--" "-." ":"];
    LineWidths = [1 1 1 1.5];
    hold on
    for k=1:length(A)
        plot(cos(deg2rad(I)),r_LSO(:,k),'Color','k','LineWidth',LineWidths(k),'LineStyle',LineStyles(k),'DisplayName',sprintf('$a = %.1f$',A(k)))
    end
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)



    %{
    if strcmp(type,'LSSO')    
        text(5.2, 0.95, '\texttt{0.0M}','Interpreter', 'latex', 'FontSize', 14); 
        text(4.45, 0.75, '\texttt{0.3M}','Interpreter', 'latex', 'FontSize', 14);
        text(3.8,0.6 ,'\texttt{0.6M}','Interpreter', 'latex', 'FontSize', 14); 
        text(2.4,0.41 ,'\texttt{a = 0.9M}','Interpreter', 'latex', 'FontSize', 14);
    else
        text(2.73,0.96, '\texttt{0.0M}','Interpreter', 'latex', 'FontSize', 14); 
        text(2.45,0.71, '\texttt{0.3M}','Interpreter', 'latex', 'FontSize', 14);
        text(2.13,0.56,'\texttt{0.6M}','Interpreter', 'latex', 'FontSize', 14); 
        text(1.6,0.29 ,'\texttt{a = 0.9M}','Interpreter', 'latex', 'FontSize', 14);
    end
    %}
return

