%------------------------------------------------------------------------------------------
% Project   :  Enhanced Modal Identification for Long-term Integrity Assessment (EMILIA)
% Author    :  Manuel F. Hormazabal
% Contact   :  m.hormazabal@gmail.com
%------------------------------------------------------------------------------------------
function EMILIA_3D_Shape_Animation_v1(Displa_Data_RAW,Bayes_PD_fn,Mode,IniSample,...
                                      EndSample,Geo_File,Make_Video,Video_Name)
%------------------------------------------------------------------------------------------
    Geo_Data  = readtable(Geo_File);
    x_coord   = Geo_Data.x_coord;
    y_coord   = Geo_Data.y_coord;
    z_coord   = Geo_Data.z_coord;
    x_chann   = Geo_Data.x_channel;
    y_chann   = Geo_Data.y_channel;
    z_chann   = Geo_Data.z_channel;
    if isempty(IniSample)
        IniSample = 1;
    end
    if isempty(EndSample)
        EndSample = numel(Displa_Data_RAW(1,:,1));
    end
    Norm_Data = Displa_Data_RAW(:,[IniSample:EndSample],Mode)./max(max(abs(Displa_Data_RAW(:,[IniSample:EndSample],Mode))));
    % ---
    
    % ---
    Phi_X  = x_coord+Norm_Data(1:numel(x_chann),1);
    Phi_Y  = y_coord+Norm_Data((numel(x_chann)+1):(numel(x_chann)+numel(y_chann)),1);  
    if isnan(z_chann(1))
        Phi_Z  = z_coord+zeros(size(z_chann));
    else
        Phi_Z  = z_coord+Norm_Data(z_chann,1); 
    end
    Lo_Marg = 2; Hi_Marg = 2;
    FigPHI = figure; FigPHI.Position = [1038,48,849,905];
    Plot_1 = scatter3(Phi_X,Phi_Y,Phi_Z,25,'filled');
    xlim([min(x_coord)-Lo_Marg max(x_coord)+Hi_Marg]);
    ylim([min(y_coord)-Lo_Marg max(y_coord)+Hi_Marg]);
    zlim([0 max(z_coord)+2]); %zlim([min(z_coord) max(z_coord)+2]);
    xlabel('Normalized Ux'); ylabel('Normalized Uy'); zlabel('Normalized Uz'); 
    grid on; grid minor; 
    title({['Mode shape Nº',num2str(Mode),'   |   MEV = ',num2str(round((Bayes_PD_fn(1,Mode)),3)),'[Hz].'],['Time step :   ',num2str(IniSample)]},'FontSize',12);
    view(45,10)
    for i = 2:numel(Norm_Data(1,:))
        if ~ishghandle(FigPHI)
            break 
        end
        Phi_X  = x_coord+Norm_Data(1:numel(x_chann),i);
        Phi_Y  = y_coord+Norm_Data((numel(x_chann)+1):(numel(x_chann)+numel(y_chann)),i);    
        if isnan(z_chann(1))
            Phi_Z  = z_coord+zeros(size(z_chann));
        else
            Phi_Z  = z_coord+Norm_Data((numel(x_chann)+2+numel(y_chann)):(numel(x_chann)+2+numel(y_chann))+numel(z_chann),i); 
        end
        set(Plot_1,'xdata',Phi_X);
        set(Plot_1,'ydata',Phi_Y);
        title({['Mode shape Nº',num2str(Mode),'   |   MEV = ',num2str(round((Bayes_PD_fn(1,Mode)),3)),'[Hz].'],['Time step :   ',num2str(IniSample+i)]},'FontSize',12);
        drawnow;
        if isequal(Make_Video,'Yes')
            Video_frames(i-1) = getframe(gcf);
        end
    end  
    if isequal(Make_Video,'Yes')
        writerObj = VideoWriter(Video_Name);
        writerObj.FrameRate = 30;
        open(writerObj);
        writeVideo(writerObj,Video_frames);
        close(writerObj); 
    end       
end
%------------------------------------------------------------------------------------------
