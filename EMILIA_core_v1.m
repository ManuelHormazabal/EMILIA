%------------------------------------------------------------------------------------------
% Project   :  Enhanced Modal Identification for Long-term Integrity Assessment (EMILIA)
% Author    :  Manuel F. Hormazabal
% Contact   :  m.hormazabal@gmail.com
%------------------------------------------------------------------------------------------
function [Acc_Data,Displa_Data_RAW,iFr,iAm,fs,Level,Channels,Bayes_PD,PDF_Fun,f_vector,...
          MEV_x,Test_Array,Output_table] = EMILIA_core_v1(Data_File,Geo_File,fs,Level,...
                                                          Decimate,Factor,TimeAlign,...
                                                          Wavelet,Distri,Dist_Res,...
                                                          Max_Damp,Min_Prob,Include_Z,...
                                                          FiltHi_Acc,FiltLo_Acc,Correct_Acc)
%------------------------------------------------------------------------------------------
%
%
%
%
%
% --- STAGE A: Geometry data import -------------------------------------------------------
    Geo_Data  = readtable(Geo_File);
    if isequal(Include_Z,'Yes')
        Channels = (cat(1,Geo_Data.x_channel,Geo_Data.y_channel,Geo_Data.z_channel))';
    else
        Channels = (cat(1,Geo_Data.x_channel,Geo_Data.y_channel))';
    end
    if isequal(Decimate,'Yes')
        fs = fs/Factor;
    end
    fprintf('\n');
    disp('  ANALISYS STARTED'); t_total = tic;
    fprintf('\n');
    disp(['    Accleration data  : ',Data_File]);               
    disp(['    Geometry File     : ',Geo_File]);
    disp(['    Selected Channels : ',num2str(Channels)]);
    disp(['    Sampling frequency: ',num2str(fs),' Hertz']);
    disp(' ');
%
% --- STAGE A: MODWPT Decomposition -------------------------------------------------------
    tic; TempT  = fprintf('    STAGE A: Computing data decomposition...');
    Raw_Data    = readmatrix(Data_File); 
    Acc_Data    = Raw_Data(:,Channels); 
    % --- Low Pass Filter
    if isequal(FiltHi_Acc,'Yes')
        fc_Hi            = 0.9*(fs/2);
        order_Hi         = 6; 
        [z_Hi,p_Hi,k_Hi] = butter(order_Hi,fc_Hi/(fs/2),'low');
        [sos_Hi,g_Hi]    = zp2sos(z_Hi,p_Hi,k_Hi); %fvtool(sos_Hi);
        for k = 1:numel(Channels)
            Acc_Data_TEMP = filtfilt(sos_Hi,g_Hi,Acc_Data(:,k));
            Acc_Data(:,k) = Acc_Data_TEMP;
        end
        clear fc_Hi order_Hi z_Hi p_Hi k_Hi sos_Hi g_Hi Acc_Data_TEMP
    end
    % --- High Pass Filter
    if isequal(FiltLo_Acc,'Yes')
        fc_Lo            = 0.1;
        order_Lo         = 6; 
        [z_Lo,p_Lo,k_Lo] = butter(order_Lo,fc_Lo/(fs/2),'high');
        [sos_Lo,g_Lo]    = zp2sos(z_Lo,p_Lo,k_Lo); %fvtool(sos_Lo);
        for k = 1:numel(Channels)
            Acc_Data_TEMP = filtfilt(sos_Lo,g_Lo,Acc_Data(:,k));
            Acc_Data(:,k) = Acc_Data_TEMP;
        end
        clear fc_Lo order_Lo z_Lo p_Lo k_Lo sos_Lo g_Lo Acc_Data_TEMP
    end
    % --- Base-line correction
    if isequal(Correct_Acc,'Yes')
        t = (0:1/fs:(numel(Acc_Data(:,1))/fs)-(1/fs))';
        fcut = (fs/2)-1;      
        for k = 1:numel(Channels)
            Vel_Data      = cumtrapz(Acc_Data(:,k))*t(2);       
            Vel_Fit       = polyfit(t,Vel_Data,2);     
            ptd           = Vel_Fit(1)*t+Vel_Fit(2);
            udd_cor1      = Acc_Data(:,k) - ptd;       
            udd_cor1_f    = fft(udd_cor1);
            df            = 1/t(2)/length(t);
            f             = df:df:df*length(t);
            [b,a]         = butter(4,fcut/(f(end)/2),'high');
            udd_cor2_f    = filter(b,a,udd_cor1_f);
            udd_cor2      = ifft(udd_cor2_f);      
            Acc_Data_Temp = real(udd_cor2);   
            Acc_Data(:,k) = Acc_Data_Temp;
        end
        clear a b df disp f fcut k ptd t u_cor ud_cor udd_cor1 udd_cor1_f udd_cor2 udd_cor2_f vel vel_Data Vel_Fit Acc_Data_Temp Vel_Data
    end
    % ---
    if isequal(Decimate,'Yes')
        Acc_Data_temp = zeros(numel(Acc_Data(:,1))/Factor,numel(Channels));
        for k = 1:numel(Channels)
            Acc_Data_temp(:,k) = decimate(Acc_Data(:,k),Factor);
        end
        Acc_Data = Acc_Data_temp;
    end
    Decomp_Data = zeros(numel(Channels),numel(Acc_Data(:,1)),2^Level);
    wt_c_freq   = zeros(2^Level,numel(Channels));
    for k = 1:numel(Channels)
        if isequal(TimeAlign,'Yes')
            [Decomp_Temp,~,wt_c_freq(:,k)] = modwpt(Acc_Data(:,k),Wavelet,Level,'TimeAlign',true); % [OUTPUT]  MODWPT data decomposition [60000 x 2^Level]
        else
            [Decomp_Temp,~,wt_c_freq(:,k)] = modwptdetails(Acc_Data(:,k),Wavelet,Level); % [OUTPUT]  MODWPT data decomposition [60000 x 2^Level]
        end
        Decomp_Data(k,:,:) = Decomp_Temp';                             
    end
    clear Decomp_Temp;
    fprintf(repmat('\b',1,TempT));
    fprintf('    STAGE A: Computing data decomposition...           DONE');
    fprintf('\n');
    disp(['           Elapsed time: ',num2str(toc),' seconds']);
    fprintf('\n');
%
% --- STAGE A: Hilbert Transform TFA ------------------------------------------------------
    tic; TempT = fprintf('    STAGE A: Computing Hilbert TFA...');
    HTZ = zeros(numel(Channels),numel(Acc_Data(:,1)),2^Level);
    iAm = zeros(numel(Channels),numel(Acc_Data(:,1)),2^Level);
    iFr = zeros(numel(Channels),numel(Acc_Data(:,1))-1,2^Level);
    for j = 1:numel(Channels)
        for k = 1:2^Level
            HTZ(j,:,k) = hilbert(Decomp_Data(j,:,k));
            iFr(j,:,k) = fs/(2*pi)*diff(unwrap(angle(HTZ(j,:,k)))); % [OUTPUT] Instantanepus Frequency [59999 x 2^Level]
            iAm(j,:,k) = abs(HTZ(j,:,k)); % [OUTPUT] Instantaneous Amplitude [60000 x 2^Level]
        end
    end
    fprintf(repmat('\b',1,TempT));
    fprintf('    STAGE A: Computing Hilbert TFA...                  DONE');
    fprintf('\n');
    disp(['           Elapsed time: ',num2str(toc),' seconds']);
    fprintf('\n');
%
% --- STAGE B: Probabilities Analysis -----------------------------------------------------
    tic; TempT = fprintf('    STAGE B: Computing Probability analysis...');
    f_vector   = 0:(fs/2)/Dist_Res:(fs/2)-((fs/2)/Dist_Res);
    PDF_Fun    = zeros(numel(Channels),Dist_Res,2^Level);
    PDF_Pb     = zeros(numel(Channels),2^Level);
    PDF_fn     = zeros(numel(Channels),2^Level);
    PDF_fn_avg = zeros(1,2^Level);
    PDF_fn_var = zeros(1,2^Level);
    for j = 1:numel(Channels)
        for k = 1:2^Level
            pdist_fit         = fitdist(iFr(j,:,k)',Distri);
            PDF_Fun(j,:,k)    = pdf(pdist_fit,f_vector); % [OUTPUT] Probability Density Functions (PDF) [(((fs/2)/(1/Res))+1) x 2^Level].
            [PDF_Pb(j,k),Loc] = max(PDF_Fun(j,:,k)); % [OUTPUT] Higher probability [1 x 2^Level].
            PDF_fn(j,k)       = f_vector(Loc); % [OUTPUT] PDF higher probability in Hertz [Channels x 2^Level].
        end  
    end
    for k = 1:2^Level
        PDF_fn_avg(k) = mean(PDF_fn(:,k)); % [OUTPUT] Averaged fn througth all channels [1 x 2^Level].
        PDF_fn_var(k) = var(PDF_fn(:,k));
    end
    fprintf(repmat('\b',1,TempT));
    fprintf('    STAGE B: Computing Probability analysis...         DONE');
    fprintf('\n');
    disp(['           Elapsed time: ',num2str(toc),' seconds']);
    fprintf('\n');
%
% --- STAGE B: Bayes Conditional Probabilities --------------------------------------------
    tic; TempT = fprintf('    STAGE B: Computing Bayes analysis...');
    Bayes_PD = zeros(Dist_Res,2^Level);
    MEV_x    = zeros(1,2^Level); 
    Var_x    = zeros(1,2^Level);
    Likehood = zeros(1,2^Level);
    for k = 1:2^Level                                                                      
        Likehood = PDF_Fun(1,:,k);   
        for i = 2:1:size(PDF_Fun,1)
            p_cond_i_1 = Likehood;                                  
            p_cond_i   = trapz(PDF_Fun(i,:,k).*p_cond_i_1)*(f_vector(2)-f_vector(1));          
            Likehood   = PDF_Fun(i,:,k).*p_cond_i_1/p_cond_i;
        end
        Bayes_PD(:,k)  = Likehood;  % [OUTPUT] Bayes likehood PDFs [2^Level x Res].
        MEV_x(k)       = trapz(f_vector.*Bayes_PD(:,k)')*(f_vector(2)-f_vector(1));                                    
        Var_x(k)       = trapz((f_vector-MEV_x(k)).^2.*Bayes_PD(:,k)')*(f_vector(2)-f_vector(1));    
    end
    fprintf(repmat('\b',1,TempT));
    fprintf('    STAGE B: Computing Bayes analysis...               DONE');
    fprintf('\n');
    disp(['           Elapsed time: ',num2str(toc),' seconds']);
    fprintf('\n');
%
% --- STAGE B: Minimum Probability & Maximum Damping TEST ---------------------------------
    if isempty(Min_Prob)
        if (Level > 2) && (Level <= 9)
            Min_Prob = 0.1*(Level+1); 
        elseif Level > 9
            Min_Prob = 0.90;
        else
            Min_Prob = 0.20;
        end
    end
    if isempty(Max_Damp)
        Max_Damp = 0.05;
    end
    Bayes_PD_Max     = zeros(1,2^Level);
    Bayes_PD_Max_Loc = zeros(1,2^Level);
    Bayes_PD_fn      = zeros(1,2^Level);
    Bayes_PD_Dp      = zeros(1,2^Level);
    for k = 1:2^Level
        [Bayes_PD_Max(k),Bayes_PD_Max_Loc(k)] = max(Bayes_PD(:,k)); 
        Half_amp = Bayes_PD_Max(k)/2;
        HBW_sel  = [];
        for j = 1:Dist_Res
            Testing_Sample = Bayes_PD(j,k);
            if (Testing_Sample >= Half_amp)
                HBW_sel = cat(1,HBW_sel,f_vector(j));
            end
        end
        Bayes_PD_fn(k) = f_vector(Bayes_PD_Max_Loc(k));        
        Bayes_PD_Dp(k) = (HBW_sel(end)-HBW_sel(1))/(2*Bayes_PD_fn(k)); 
    end
    Bayes_PD_Norm = Bayes_PD./max(max(abs(Bayes_PD)));
    Test_Array = ones(2^Level,1);
    for k = 1:2^Level
        if (max(Bayes_PD_Norm(:,k)) >= Min_Prob)...
        && (Bayes_PD_Dp(k) <= Max_Damp)
            Test_Array(k) = 2; % [OUTPUT] Test results for all components [1 x 2^Level].   
        end
    end
%    
% --- STAGE B: Compute Displacements ------------------------------------------------------
    tic; TempT = fprintf('    STAGE B: Computing displacements...');
    t0              = (0:1/fs:(numel(Acc_Data(:,1))/fs)-(1/fs))';
    [order_Hi,Wc_Hi]      = buttord(0.0001,0.9,0.001,60);
    [bHi,aHi]         = butter(order_Hi,Wc_Hi,'high');
    Displa_Data_RAW = zeros(size(Decomp_Data));
    for k = 1:2^Level
        for j = 1:numel(Channels)
            Mode_Vel               = cumtrapz(t0,Decomp_Data(j,:,k));       
            Filt_Vel               = filtfilt(bHi,aHi,Mode_Vel);
            Non_Filt_Displa        = cumtrapz(t0,Filt_Vel);
            Displa_Data_RAW(j,:,k) = filtfilt(bHi,aHi,Non_Filt_Displa); % [OUTPUT] Mode Shapes: Displacement Hi-Pass filtering [Channels,Samples,2^Level]            
        end
    end
    Displa_Data_NORM = Displa_Data_RAW/...
                       max(max(max(abs(Displa_Data_RAW)))); % [OUTPUT] Mode Shapes: Normalized Displacement [Channels,Samples,2^Level] 
    fprintf(repmat('\b',1,TempT));
    fprintf('    STAGE B: Computing displacements...                DONE');
    fprintf('\n');
    disp(['           Elapsed time: ',num2str(toc),' seconds']);
    fprintf('\n');   
%
% --- STAGE B: Analysis Output ------------------------------------------------------------
    Output_table = table((1:2^Level)',round(PDF_fn_avg,2)',...
           round(PDF_fn_var,3)',round(MEV_x,2)',round(Var_x,3)',round(Bayes_PD_Dp,3)',Test_Array);
    Output_table.Properties.VariableNames = {'Nº','MEV(Avg)','Var(Avg)','MEV(Bayes)','Var(Bayes)','Damping(Bayes)','Test'}; 
    fprintf('\n');
    disp('    ANALISYS COMPLETED'); 
    disp(['      Total elapsed time: ',num2str(toc(t_total))]);
    fprintf('\n');
    fprintf('\n');
    disp(' ');
    disp('   ------------------------------- ANALYSIS RESULTS --------------------------------');
    disp(' ');
    disp(Output_table);
    disp('   ---------------------------------------------------------------------------------');
    disp(' ');
    fprintf('\n');
    fprintf('\n');
    disp('  ANALISYS COMPLETED'); t_total = tic;
    fprintf('\n');
    fprintf('\n');
    fprintf('\n');
    fprintf('\n');
%
