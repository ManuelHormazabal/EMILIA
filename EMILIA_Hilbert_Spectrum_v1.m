%------------------------------------------------------------------------------------------
% Project   :  Enhanced Modal Identification for Long-term Integrity Assessment (EMILIA)
% Author    :  Manuel F. Hormazabal
% Contact   :  m.hormazabal@gmail.com
%------------------------------------------------------------------------------------------
function EMILIA_Hilbert_Spectrum_v1(Acc_Data,iFr,iAm,fs,Level,Channels,Ch_plot,Sq_plot,...
                                    Sp_skip,freq_span,time_span,Theme,Order,Res,Wavelet,...
                                    Mov_AVG,Mov_WIN,Mov_skip)
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
    if isequal(Theme,'Blue') || isempty(Theme)
        Background = [0 0 0.53];
        Mov_Color  = 'w';
    end
    if isequal(Theme,'White')
        Background = [1 1 1];
        Mov_Color  = 'k';
    end
    t0 = 0:1/fs:(numel(Data_iAm(1,:,1))/fs)-(1/fs);
    if isequal(Order,'Global')
        [~,Loc] = sort((rms(Acc_Data)));
    end
    for i = 1:numel(Sq_plot)
        if isequal(Order,'Local') 
            [~,Loc] = sort(rms(Data_iAm(:,:,i)'));
        end
        for k = 1:numel(Ch_plot)   
            hold on; 
            scatter(t0(2:Sp_skip:end),Data_iFr(Loc(k),1:Sp_skip:end,Sq_plot(i)),Res,Data_iAm(Loc(k),2:Sp_skip:end,Sq_plot(i)),'filled'); 
        end
        if isequal(Mov_AVG,'Yes')
            line(t0(2:Mov_skip:end),movmean(Data_iFr(Loc(k),1:Mov_skip:end,Sq_plot(i)),Mov_WIN),'linewidth',1,'Color',Mov_Color);
        end
    end
    xlabel('Time [s]'); 
    ylabel('Instantaneous Frequency [Hz]');       
    title({['Hilbert Spectrum - Wavelet: ',Wavelet],[num2str(Level),'th level MODWPT decomposition']}); 
    barra_T = colorbar; barra_T.FontSize = Font_Size; 
    barra_T.Label.String = 'Instantaneous Amplitude'; barra_T.Label.FontSize = Font_Size; 
    colormap jet(128);
    if isempty(time_span)
        time_span  = [0 t0(end)]; 
    end
    if isempty(freq_span)
        freq_span  = [0 fs/2]; 
    end
    ax = gca; 
         ax.XLim     = time_span;
         ax.YLim     = freq_span; 
         ax.Color    = Background;  
         ax.FontSize = Font_Size; 
         
end
%------------------------------------------------------------------------------------------
