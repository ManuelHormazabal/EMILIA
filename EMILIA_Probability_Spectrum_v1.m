%------------------------------------------------------------------------------------------
% Project   :  Enhanced Modal Identification for Long-term Integrity Assessment (EMILIA)
% Author    :  Manuel F. Hormazabal
% Contact   :  m.hormazabal@gmail.com
%------------------------------------------------------------------------------------------
function EMILIA_Probability_Spectrum_v1(Ch_plot,Sq_plot,Sp_skip,iFr,iAm,Bayes_PD,...
                                        PDF_Fun,Wavelet,Level,fs,Channels,f_vector,...
                                        MEV_x,Test_Array,Min_Prob,Y_Axis,Plot_PDs,...
                                        Plot_BPD,Plot_Test,fn_span)
%------------------------------------------------------------------------------------------
    Font_Size = 18;
    if isempty(Sq_plot)
        Sq_plot = 1:2^Level;
    end
    if isempty(Ch_plot)    
        Ch_plot = 1:numel(Channels);
    end
    Data_iFr = iFr(Ch_plot,:,:);
    Data_iAm = iAm(Ch_plot,:,:);
    figure('units','normalized','outerposition',[0 0 1 1]);
    tt = 0:(1/fs):(length(Data_iAm(1,:,1))/fs)-(1/fs);
    if isempty(fn_span)
        fn_span = [0 fs/2];
    end
    xlim(fn_span); 
    grid on; grid minor; xlabel('Instantaneous Frequency [Hz]');
    ax = gca; 
        ax.YColor = 'r'; 
        ax.FontSize = Font_Size;
    if isempty(Sp_skip) == 0
        c = colorbar('northoutside','FontSize',Font_Size);
        c.Label.String = 'Time [s]'; c.Label.FontSize = Font_Size; 
        colormap parula(128); 
    end
    yyaxis left
        if isempty(Sp_skip) == 0
            for j = 1:numel(Ch_plot)
                for k = 1:numel(Sq_plot)
                    hold on;
                    S1 = scatter(squeeze(Data_iFr(j,1:Sp_skip:end,Sq_plot(k))),squeeze(Data_iAm(j,2:Sp_skip:end,Sq_plot(k))),1,tt(2:Sp_skip:end),'filled','Displayname','Data Samples');
                end
            end
        end
        ylabel('Instantaneous Amplitude');
    yyaxis right
        L1 = []; L2 = [];
        if isequal(Plot_PDs,'Yes')
            for j = 1:numel(Ch_plot)
                for k = 1:numel(Sq_plot)
                    if k == 1 && j == 1
                        L1 = line(f_vector,PDF_Fun(j,:,Sq_plot(k)),'Color',[0.6350 0.0780 0.1840],'LineStyle','--','LineWidth',0.8,'Displayname','Kernel - Avg. fn = 1.95 Hz.');
                    else
                        line(f_vector,PDF_Fun(j,:,Sq_plot(k)),'Color',[0.6350 0.0780 0.1840],'LineStyle','--','LineWidth',0.8,'HandleVisibility','off');
                    end
                end
            end
        end
        if isequal(Plot_BPD,'Yes')
            for k = 1:numel(Sq_plot)
                if k == 1
                    L2 = line(f_vector,Bayes_PD(:,Sq_plot(k)),'Color','r','LineStyle','-','LineWidth',3,'Displayname','Bayes - MEV fn = 1.87 Hz.');
                else
                    line(f_vector,Bayes_PD(:,Sq_plot(k)),'Color','r','LineStyle','-','LineWidth',3,'HandleVisibility','off');
                end
            end
        end
        ylabel('Probability Density');
        if isequal(Plot_Test,'Yes')
            Lbl_cont = 1;
            for k = 1:numel(Sq_plot) 
                if Test_Array(Sq_plot(k)) == 2
                    Lbl_cont = Lbl_cont + 1;
                    LB(Lbl_cont) = line(f_vector,Bayes_PD(:,Sq_plot(k)),'Color',[0.4 1-((k/2)/(2^Level)) 0.2],'LineStyle','--','LineWidth',2);
                    T1{Lbl_cont} = ['MEV Nº',num2str(Lbl_cont-1),':   ',num2str(round(MEV_x(Sq_plot(k)),3)),' Hz.'];
                end   
            end
            LB(1) = line([0 fs/2],[Min_Prob*(max(max(abs(Bayes_PD)))) Min_Prob*(max(max(abs(Bayes_PD))))],'Color','b','LineStyle','--','LineWidth',0.8);
            T1{1} = ['Prob. Threshold (',num2str(100*Min_Prob),'%)'];
            leg2 = legend(LB,T1,'Location','NorthEast','FontSize',Font_Size); 
            title(leg2,{['TEST RESULTS']}); 
        end
        if isequal(Y_Axis,'Log')
            yyaxis left
                set(gca,'YScale','log'); 
                ylim([(max(max(max(iAm)))/200) (20*max(max(max(iAm))))]);            
            yyaxis right
                set(gca,'YScale','log'); 
                ylim([(max(max(Bayes_PD))/200) (20*max(max(Bayes_PD)))]);
        end
        if isempty(Sp_skip)
            title({['EMILIA Probability Spectrum'],[num2str(Level),'th level ',Wavelet,' MODWPT decomposition'],['Input Channels:  ',num2str(Channels)]});
        else     
            title({['EMILIA Time-Dependent Probability Spectrum'],[num2str(Level),'th level ',Wavelet,' MODWPT decomposition'],['Input Channels:  ',num2str(Channels)]});
        end
        if isequal(Plot_Test,'Yes')
            aleg = axes('position',get(gca,'position'),'visible','off');
        else
            aleg = gca;
        end
        if isempty(Sp_skip) == 0 && isempty(L1) == 0 && isempty(L2) == 0
            leg1 = legend(aleg,[S1,L1,L2],'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');  
        elseif  isempty(Sp_skip) && isempty(L1) == 0 && isempty(L2) == 0
            leg1 = legend(aleg,[L1,L2],'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');
        elseif  isempty(Sp_skip) && isempty(L1) == 0 && isempty(L2)
            leg1 = legend(aleg,L1,'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');
        elseif  isempty(Sp_skip) && isempty(L1) && isempty(L2) == 0
            leg1 = legend(aleg,L2,'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');  
        elseif isempty(Sp_skip) == 0 && isempty(L1) == 0 && isempty(L2)
            leg1 = legend(aleg,[S1,L1],'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');
        elseif isempty(Sp_skip) == 0 && isempty(L1) && isempty(L2) == 0
            leg1 = legend(aleg,[S1,L2],'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');
        elseif isempty(Sp_skip) == 0 && isempty(L1) && isempty(L2)
            leg1 = legend(aleg,S1,'Location','NorthWest','FontSize',Font_Size);
            title(leg1,'PDFs');
        end
end
%------------------------------------------------------------------------------------------
