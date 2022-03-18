%------------------------------------------------------------------------------------------
% Project   :  Enhanced Modal Identification for Long-term Integrity Assessment (EMILIA)
% Author    :  Manuel F. Hormazabal
% Contact   :  m.hormazabal@gmail.com
%------------------------------------------------------------------------------------------
function EMILIA_2D_Shape_Animation_v1(Displa_Data_RAW,Bayes_PD_fn,Mode,IniSample,...
                                      EndSample,Geo_File,Make_Video,Video_Name)
%------------------------------------------------------------------------------------------
    Geo_Data  = readtable(Geo_File);
    x_coord   = Geo_Data.x_coord;
    y_coord   = Geo_Data.y_coord;
    x_chann   = Geo_Data.x_channel;
    y_chann   = Geo_Data.y_channel;
    if isempty(IniSample)
        IniSample = 1;
    end
    if isempty(EndSample)
        EndSample = numel(Displa_Data_RAW(1,:,1));
    end
    Norm_Data = Displa_Data_RAW(:,[IniSample:EndSample],Mode)./max(max(abs(Displa_Data_RAW(:,[IniSample:EndSample],Mode))));
    FigPHI = figure; FigPHI.Position = [570,418,1304,535];
    Phi_X  = Norm_Data(x_chann,1);
    Phi_Y  = Norm_Data(y_chann,1);  
    Plot_1 = line(x_coord+Phi_X,y_coord+Phi_Y,'linewidth',2); 
    hold on; 
    Plot_2 = scatter(x_coord,y_coord,25,'filled');
    xlim([min(x_coord)-1 max(x_coord)+1]);
    ylim([-1 1]);
    xlabel('Normalized Ux'); ylabel('Normalized Uy'); 
    grid on; grid minor; 
    title({['Mode shape Nº',num2str(Mode),'   |   MEV = ',num2str(round((Bayes_PD_fn(1,Mode)),3)),'[s].'],['Time step :   ',num2str(IniSample)]},'FontSize',12); 
    for i = 2:numel(Norm_Data(1,:))
        if ~ishghandle(FigPHI)
            break
        end
        Phi_X  = x_coord+Norm_Data(x_chann,i);
        Phi_Y  = Norm_Data(y_chann,i);  
        set(Plot_1,'xdata',Phi_X);
        set(Plot_1,'ydata',Phi_Y);
        set(Plot_2,'xdata',Phi_X);
        set(Plot_2,'ydata',Phi_Y);
        title({['Mode shape Nº',num2str(Mode),'   |   MEV = ',num2str(round((Bayes_PD_fn(1,Mode)),3)),' Hz.'],['Time-Step :   ',num2str(IniSample+i)]},'FontSize',12);
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
