function twisted_wf(l,m,aeq,source)
    % ========================================================================================================================
    % Produces Figs. 17-18 of Paper I
    % (l,m) = multipolar indices of waveform mode
    % aeq = BH spin for the equatorial waveform
    % source = can be either 'teukode' or 'data'
    % ========================================================================================================================

    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    norm = sqrt((l+2).*(l+1).*l.*(l-1));

    if strcmp(source,'teukode')
        % Reads waveform modes from Teukode output
        % Rotates waveform, writes to file the equatorial rotated and the precessing waveforms in the py/ directory 
        % To align run the script found in the py directory, then copy the outputs in data/ and run this again with source = 'data'

        write = 1;

        if aeq==0.35
            load('~/waveforms/K/twist/i0/a035/wf.mat')
            seq = s;
            load('~/waveforms/K/twist/i45/wf.mat')
            sp = s;
        elseif aeq==0.2
            load('~/waveforms/K/twist/i0/a02/wf.mat')
            seq = s;
            load('~/waveforms/K/twist/i30/a023/wf.mat')
            sp = s;
        else 
            error('no configurations to twist with the given spin')
        end


        sr = struct;
        sr = DB_mode_rotate_timedep(seq,sp.dyn,l,m,sr,'direction','backward');
        %sr = DB_mode_rotate_upgrade(seq,sp.dyn,l,m,sr);

        tLR_eq = tLR_splined(seq.dyn);
        tLR_prec = tLR_splined(sp.dyn);
        dt = tLR_eq - tLR_prec;
        if m>=0
            t = sp.ell(l).emm(m+1).t;
            hlm_prec = sp.ell(l).emm(m+1).hlm./norm;
            hlm_rot = spline(sr.ell(l).emm(m+1).t - dt,sr.ell(l).emm(m+1).hlm./norm,t);
        else
            t = sp.ell(l).emminus(-m+1).t;
            hlm_prec = sp.ell(l).emminus(-m+1).hlm./norm;
            hlm_rot = spline(sr.ell(l).emminus(-m+1).t - dt,sr.ell(l).emminus(-m+1).hlm./norm,t);
        end

        t_in = tLR_prec - 500;
        t_end = tLR_prec + 50;

        figure
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,abs(hlm_prec),'LineWidth',my_linewidth,'LineStyle','-','Color','b','DisplayName','$|\Psi_{22}^{\rm prec}|/\nu$')
        plot(t,abs(hlm_rot),'LineWidth',my_linewidth,'LineStyle','-','Color','r','DisplayName','$|\tilde{\Psi}_{22}^{\rm eq}|/\nu$')
        plot(t,real(hlm_prec),'LineWidth',my_linewidth,'LineStyle','-','Color',[0 0 1 0.3],'DisplayName','$\Re[\Psi_{22}^{\rm prec}]/\nu$')
        plot(t,real(hlm_rot),'LineWidth',my_linewidth,'LineStyle','-','Color',[1 0 0 0.3],'DisplayName','$\Re[\tilde{\Psi}_{22}^{\rm eq}]/\nu$')
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
        xlim([t_in,tLR_prec])
        lims = ylim;
        xlim([t_in,t_end])
        ylim(lims);

        basedir = '~/repos/teobiresumsprecessing/py/';
        rot_name = sprintf('aeq_%.2f_lm_%d%d_rotated.dat',aeq,l,m);
        prec_name = sprintf('aeq_%.2f_lm_%d%d_precessing.dat',aeq,l,m);

        if write
            fprintf('Writing to %s ... ',basedir)
            fid1 = fopen([basedir,rot_name],'w');
            fid2 = fopen([basedir,prec_name],'w');

            for i=1:length(t)
                fprintf(fid1,'%15.12f %15.12f %15.12f \n',t(i),real(hlm_rot(i)),imag(hlm_rot(i)));
                fprintf(fid2,'%15.12f %15.12f %15.12f \n',t(i),real(hlm_prec(i)),imag(hlm_prec(i)));
            end

            fclose(fid1);
            fclose(fid2);
            fprintf('Done.\n')
        end

    elseif strcmp(source,'data')
        % read in waveforms (rotated equatorial and inclined) aligned via python script
        % output of the script must be placed in Paper/data

        if aeq==0.35
            delta=1;
        else
            delta=0;
        end
        basedir = 'Paper/data/';
        rot_name = sprintf('aeq_%.2f_lm_%d%d_rotated_shifted.dat',aeq,l,m);
        prec_name = sprintf('aeq_%.2f_lm_%d%d_precessing_shifted.dat',aeq,l,m);

        data = readmatrix([basedir,prec_name]);
        t_p = data(:,1);
        hlm_p = data(:,2) + 1i.*data(:,3);
        clear data

        data = readmatrix([basedir,rot_name]);
        t_r = data(:,1);
        hlm_r = data(:,2) + 1i.*data(:,3);

        idx_end = find(abs(hlm_p)<1e-8,1);
        t_end = t_p(idx_end);
        t_in = t_end-1000;

        phase_p = -unwrap(angle(hlm_p));
        phase_r = spline(t_r,-unwrap(angle(hlm_r)),t_p);
        A_p = abs(hlm_p);
        A_r = spline(t_r,abs(hlm_r),t_p);
        DeltaPhi = phase_r - phase_p - 2*pi.*delta;
        DeltaA = abs(A_r-A_p)./(A_p);   

        figure
        tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');

        nexttile(tl,[2,1])
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
        xlim([t_in,t_end])
        hold on
        plot(t_p,real(hlm_p),'LineWidth',my_linewidth,'LineStyle','-','Color','k','DisplayName','$\Re[\Psi_{22}^{\rm prec}]/\nu$')
        plot(t_r,real(hlm_r),'LineWidth',my_linewidth,'LineStyle','--','Color',[1 0 0 ],'DisplayName','$\Re[\tilde{\Psi}_{22}^{\rm eq}]/\nu$')
        plot(t_p,abs(hlm_p),'LineWidth',my_linewidth,'LineStyle','-','Color',[0 0 0 0.2],'DisplayName','$|\Psi_{22}^{\rm prec}|/\nu$')
        plot(t_r,abs(hlm_r),'LineWidth',my_linewidth,'LineStyle','--','Color',[1 0 0 0.2],'DisplayName','$|\tilde{\Psi}_{22}^{\rm eq}|/\nu$')
        ylim([-0.3,0.4])

        nexttile
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
        xlim([t_in,t_end])
        hold on
        plot(t_p,DeltaA,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
        plot(t_p,DeltaPhi,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
        if aeq==0.35
            ylim([-0.2,3.2])
        else
            ylim([-0.3,1.2])
        end
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','latex')
    end
return