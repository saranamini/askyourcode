%% *Pouch Cell FBG analysis*
% elab - documentation: |[Link to wiki about sensors and general info]|
% 
% _TODO_
%% 
% * _In phaseIdentifierGen, it is hard-coded whether it assumes the cycles under 
% consideration starts on a charge or discharge cycle. This can be problematic, 
% and an algorithmic way of determining this should be used instead. This is at 
% around line 32 in the code._
% * _Within phaseIdentifierGen(), in removeOutlierGroups(), the 'minGroupSize' 
% variable can cause problems; for high C rates, it can be much smaller (~1000), 
% but for longer tests with low C-rates, it sometimes needs to be over 7000 to 
% remove incomplete phases. Once again, this should be determined algorithmically._
% * _Add Overall_Cali to this main script, and have it run depending on user 
% preference. Also, consider saving enough variables from this master script to 
% the .mat file so that Overall_Cali can run by itself after importing the desired 
% .mat file._
%% 
% _NOTES:_
%% 
% * _On 17 Aug 2023, Luke changed the definition of SOC in phaseIdentifierGen. 
% Because this may be a temporary change, the saved data using this definition 
% has been denoted with '..._All_Data_x*_REALKAP*.mat'. You can go into phaseIdentifierGen 
% to change this back, if you'd like. Change both the location of the definition 
% (line ~89), and where it saves (at the bottom of the function script). For example, 
% realKap=3.1900Ah for C/25, and is 2.8787 Ah for 1C._
% * _The timetable "relevantDataTT" has been created to more simply look at 
% and call upon the fibre and temperature data under consideration. It makes the 
% code more general and directly relates the two datasets rather than needing 
% to call upon "S3F3" and "x_cell_temp" etc._
%% *1. Raw Data Processing*
% This section includes the initial setup and preprocessing of raw data. It 
% involves setting up paths, defining test names, and preparing the environment 
% for data processing.
%% 
% * Define Test Name and Paths
% * Test name definition (|Test_name|)
% * Path definitions 
% * Preprocessing of Raw Data
% * Handling Hero data channels (|Hero_chnl|)
% * Checking for Hero data existence and processing 

% Define Test name to be analysed
Test_name = '20230828#25degC_lowCRateTest#Ko10-xx_BLB1-xB'

% Define paths
load('Path_Var.mat')
addpath(genpath(Git_path.Sara_work))
test_path = fullfile(Data_path.Sara_work, '/Measurements/', Test_name);

% Preprocessing of raw data
FiSens_path = [test_path, '\FiSens'];
Pico_path = [test_path, '\Pico'];
Pt100_path = [test_path, '\Pt100'];

Hero_chnl = ["264", "261", "263"];
Cell_id = ["BLB1-1B", "Ko10-42", "Ko10-43"];

if ~exist([test_path, '\Hero'], 'dir')
    disp('no hero data found')
else
    path_py = strrep(test_path, '\', '/');
    for i = 1:length(Hero_chnl)
        Hero_path{i} = [path_py, '/Hero/Hero_', char(Cell_id(i))];
        if exist([test_path, '\Hero\Hero_', char(Cell_id(i)), '\Hero_complete.csv'], "file") == 2
            disp([test_path, '\Hero\Hero_', char(Cell_id(i)), '\Hero_complete.csv exists'])
            continue
        else
            pyrunfile("C:\Users\saja10\Documents\GitHub\code_overview\Python\Tools\Hero_CSV_cat_fnc_new.py", "" + ...
                "z", pth=Hero_path{i}, chl=char(Hero_chnl(i)));
        end
    end
end

%% *2. Data Import*
% This section deals with importing various data sets into the Matlab environment. 
% It includes importing data from different sources and formats.
%% 
% * FiSens data import and processing
% * Pico data import and processing 
% * Hero data import and synchronization
% * Putty data import and processing 

% FiSens Import
FiSens_S_Dis = FiSens_S_Dis_Identifier(FiSens_path);
if exist("FiSens_tt", "var")
    disp('Imported FiSens TimeTable is used')
else
    if isfile(fullfile(FiSens_path, 'FiSens_complete.mat'))
        load(fullfile(FiSens_path, 'FiSens_complete.mat'))
    else
        FiSens_tt = Import_FiSens(FiSens_path, FiSens_S_Dis, 'ISO');
    end
end
%% 
% 

% Pico Import
if ~exist("Pico_TT", "var")
    if isfile(fullfile(Pico_path, '*.mat'))
        load(fullfile(Pico_path, '*.mat'));
    else
        Pico_TT = Import_pico(Pico_path, {'T_pico_mean_cell_01 / °C', 'T_pico_mean_cell_02 / °C', 'T_pico_mean_ambient / °C'});
    end
end

%% *3. Synchronization of Data Sets*
% This section focuses on synchronizing the imported data sets to ensure they 
% are aligned and can be analyzed together.
%% 
% * Synchronize FiSens and Pico Data
% * Timezone adjustments and synchronization (|FiSens_tt.Time.TimeZone|, |synchronize|)
% * Synchronize Complete Data Sets
% * Merging and synchronizing all data sets

% Synchronize FiSens and Pico Data
FiSens_tt.Time.TimeZone = Pico_TT.Time.TimeZone;
TT_sync_FiPico = synchronize(FiSens_tt, Pico_TT, 'commonrange', 'linear');

% Load in Hero Data and Synchronize
HeroID = 1;
% if Hero_data
    if ~exist("Hero_TT", "var")
        if isfile(fullfile(Hero_path{HeroID}, 'Hero_complete.mat'))
            load(fullfile(Hero_path{HeroID}, 'Hero_complete.mat'))
        else
            Hero_TT = Import_hero(Hero_path{HeroID});
        end
    end
    Hero_TT.Time.TimeZone = Pico_TT.Time.TimeZone;
    TT_sync_com = synchronize(TT_sync_FiPico, Hero_TT, 'commonrange', 'linear');
% end

% Import and Synchronize pt100 Data
if exist(Pt100_path, "dir") == 7
    if ~exist("Pt100_TT", "var")
        if isfile(fullfile(Pt100_path, '*.mat'))
            load(fullfile(Pt100_path, '*.mat'));
        else
            Pt100_TT = Import_influxdb_pt100(Pt100_path, initial_chamber_temp);
        end
        Pt100_TT.Time = Pt100_TT.Time + Wago_ts;
    end
    if Hero_data
        TT_sync_im = synchronize(TT_sync_FiPico, Pt100_TT, 'commonrange', 'linear');
    else
        TT_sync_com = synchronize(TT_sync_FiPico, Pt100_TT, 'commonrange', 'linear');
    end
end

% Import Putty Data
if ~exist("Putty_TT", "var")
    Import_PuTTY(test_path, puttyPath, TT_sync_com, Test_name);
    load(fullfile(test_path + "\Pt100\", "Putty_complete" + ".mat"), 'Putty_TT')
end

completeTimetable = synchronize(TT_sync_com, Putty_TT, 'last', 'linear');

%% *4. Saving Data Sets in Measurement Folder*
% This section is dedicated to saving the processed and synchronized data into 
% specific formats and locations .
% 
% 
%% 
% * Saving complete timetable 
% * Create and Save Relevant Data Timetable
% * Extracting and saving relevant data

% Create and Save Relevant Data Timetable
datasetNames = "S3F3";
tempNames = "x_cell_Temp";
relevantDataTT = timetable(completeTimetable.Time, completeTimetable.S1F1);
for i = 1:numel(datasetNames)
    relevantDataTT.(datasetNames(i)) = completeTimetable.(datasetNames(i));
end
relevantDataTT = removevars(relevantDataTT, 'Var1');
for i = 1:numel(datasetNames)
    if tempNames(i) == "DNE"
        continue
    end
    relevantDataTT.(string(datasetNames(i) + "_T")) = completeTimetable.(tempNames(i));
end

save(fullfile(test_path, "completeTimetable" + Test_name + ".mat"), 'completeTimetable')


disp("done of completeTimetable")
%% 
% % Error handling and special case management

error('Script until here for data processing. Afterwards options for visualization. Delete this line if you like to test the following segments')

% Exceptions to catch special cases due to experimental circumstances
% if strcmp(Test_name, '25155015degC2Pyr-0C_Ko10-xx_F1S2F1S4_20221219')    
%     TT_sync_com = removevars(TT_sync_com, 19); % Dirty fix
%     TT_sync_com.Properties.VariableNames{18} = 'T_pico_mean_hall / °C'; % Dirty fix
% end
%% Visualization
% TODO
%% 
% * Add more documentation to each point
% * Put most of the contained code blocks into scrips that can be called from 
% the live script
% * The live script should be well documented collection of functions
% Load in of data

data_file_info = dir(fullfile(Data_path.Andre_work,'/Measurements/',Test_name,'*.mat'));
if ~exist('completeTimetable','var')
load(fullfile(Data_path.Andre_work,'/Measurements/',Test_name,data_file_info.name));
end
disp('Cycling data of cell x') %extract the information of cell in the raw data here - replace x with cell ID
% Identify relevant identifiers for test results
% This still depends on how the exactly the data is saved (is changed in paralell).
% 
% TODO
%% 
% * Get relevant identifiers like Cell_ID, date, num of fibers, sensors
% * If in doubt check with elab

%might take some work later here
TT_sync_com=completeTimetable; %needed right now because all function downbelow work for this tt
Hero_data=1; % might not be needed
%introduce relevant_data tt?

% General settings and zero reference
% General settings for visualization

zero_ref=600:620; % reference to 0
GCAfontsize=20; % fontsize
fig_siz=[100 100 1015 571]; % Figure Size - aim for target size

% In case that the relaxed state is not at the beginning of the measurement
Offset_wl.S1F2=-0.03; 
Offset_wl.S1F3=-0.007;
Offset_wl.S2F2=-0.023;

%% 
% Might be reavitvated later.

% Correction for fluctuaitions of ambient temperature of FiSens
% Ammount of correction is based on internal FBG of each channel
% Correction amount may differ depending on what the FBG is fixed on (e.g.
% pouch foil, aluminium substate, Carbon Fibre material)


% VarNames_com=TT_sync_com.Properties.VariableNames;
% num_sen=length(Meta.Fibres);
% Sensors_fibres=cell(1,num_sen);
% for i = 1:num_sen
%     Sensors_fibres{i}=sprintf('S%i',i);
%     TF=contains(VarNames_com,Sensors_fibres{i});
%     FiSens_ref_FBG=find(TF,1,'last');
%     Ambient_fixed_corr.(Sensors_fibres{i})=TT_sync_com.(VarNames_com{i})-mean(TT_sync_com.(VarNames_com{i})(zero_ref));
% end

% Create colors

num_col=3;
x=1/num_col:1/num_col:1;
for i=1:num_col
    Black{i}=[x(i), x(i), x(i)];
    Blue{i}=[0.1 0.1 x(i)];
    Red{i}=[x(i) 0.1 0.1];
end

%Custom Color Set
CusCol{1}='#FF9000';
CusCol{2}='#090C9B';
CusCol{3}='#D81159';
CusCol{4}='#0A0908';
CusCol{5}='#3C91E6';

%Tu Colors
TuCol{1}=[0 140 79]/255;
TuCol{2}=[128 128 128]/255;
TuCol{3}=[140 28 0]/255;

%Shaded colors
if strcmp(Test_name,'20230714_25degC_twenty1CRateTest_BLB1-xB')
    num_col=32;
    light=26;
    dark=14;
else
    num_col=5;
    light=5;
    dark=3;
end

for i = 1:num_col
    red{i}=[1/(num_col+1)*i 0 0];
    green{i}=[0 1/(num_col+1)*i 0];
    blue{i}=[0 0 1/(num_col+1)*i];
    margenta{i}=[1/(num_col+1)*i 0 1/(num_col+1)*i];
    cyan{i}=[0 1/(num_col+1)*i 1/(num_col+1)*i];
    yellow{i}=[1/(num_col+1)*i 1/(num_col+1)*i 0];
    grey{i}=[1/(num_col+1)*i 1/(num_col+1)*i 1/(num_col+1)*i];
end

colorList = {red{light}, margenta{light}, cyan{light}, blue{light}, yellow{light}, grey{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}, green{light}};

colorListDark={red{dark}, margenta{dark}, cyan{dark}, blue{dark}, yellow{dark}, grey{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}, green{dark}};

% Determine limits of plots

%might not be necessary for overview, but for sure looking at details.
%% FiSens wavelength shifts
% The FiSens data is saved as SxFy where x is the channel of the spectrometre 
% and y is the sensor number. The sensors are enumerated from the channel connection 
% to the end of the fibre. The sensors on the graphics are however sorted in the 
% very different direction. Therefore attention must be paid to the graphical 
% reference.
% 
% Usually this means here that is the naming of the sensors is opposite to the 
% saved format.
% 
% Example: S1F5 is the first sensors S1 if the fibre has overall 5 sensors on 
% it

fig=figure("Name",'FBG wavelength shift','DefaultAxesFontSize',GCAfontsize);

% gca.fontsize=GCAfontsize;
%Voltage and current during 0.5 C charging phase and break
if Hero_data==1
    tiledlayout(2,1,'TileSpacing','Compact')
else
    tiledlayout(1,1,'TileSpacing','Compact')
end

fig_handle=gcf;
set(gcf,'Visible','on')
set(gcf,'position',fig_siz)
set(gcf,'color','w');
set(0, 'DefaultLineLineWidth', 1.5);


ax1 = nexttile;
ax_handle=gca;

hold on;
grid on;
% xlabel('$\mathrm{t_{exp} \; / \; HH:MM:SS}$','Interpreter','latex')
ylabel('$\Delta\mathrm{\lambda \; / \; nm}$','Interpreter','latex')
% ylim([-0.05 0.25])

%(1:867867) only for this set for demonstration
plot(TT_sync_com.Time(1:867867),TT_sync_com.S1F6(1:867867)-mean(TT_sync_com.S1F6(zero_ref)),'DisplayName','S1_{T, \epsilon, Ko42}','Color',margenta{5});
plot(TT_sync_com.Time(1:867867),TT_sync_com.S1F5(1:867867)-mean(TT_sync_com.S1F5(zero_ref)),'DisplayName','S2_{T,Ko42}','HandleVisibility','on','Color',margenta{3});
plot(TT_sync_com.Time(1:867867),TT_sync_com.S1F4(1:867867)-mean(TT_sync_com.S1F4(zero_ref)),'DisplayName','S3_{T,Ko42}','HandleVisibility','on','Color',grey{3});
plot(TT_sync_com.Time(1:867867),TT_sync_com.S2F1(1:867867)-mean(TT_sync_com.S2F1(zero_ref)),'DisplayName','S1_{T,Ko43}','HandleVisibility','on','Color',grey{1});
% plot(TT_sync_com.Time,TT_sync_com.S1F3-mean(TT_sync_com.S1F3(zero_ref)),'DisplayName','S4_{T,ce1}','HandleVisibility','on','Color',cyan{4}); %+Offset_wl.S1F3
% plot(TT_sync_com.Time,TT_sync_com.S1F2-mean(TT_sync_com.S1F2(zero_ref)),'DisplayName','S5_{T, \epsilon, ce1}','Color',margenta{5}); %+Offset_wl.S1F2
% plot(TT_sync_com.Time(1:867867),TT_sync_com.S3F3(1:867867)-mean(TT_sync_com.S3F3(zero_ref)),'DisplayName','S3_{T, \epsilon, BLB1}','HandleVisibility','on','Color',green{3});

% plot(TT_sync_com.Time,TT_sync_com.S1F7-mean(TT_sync_com.S1F7(zero_ref)),'DisplayName','S_{ref,ch1}','Color','k');
% plot(TT_sync_com.Time,TT_sync_com.S1F5-mean(TT_sync_com.S1F5(zero_ref)),'DisplayName','S2_{T}','Color','k');
% plot(TT_sync_com.Time,ficorr.S1F6-mean(TT_sync_com.S1F6(zero_ref)),'DisplayName','S1_{\epsilon}','Color','c');


% plot(TT_sync_com.Time,TT_sync_com.S2F5-mean(TT_sync_com.S2F5(zero_ref)),'Color',CusCol{1},'DisplayName','S1(T3)_{T}');

% plot(TT_sync_com.Time,TT_sync_com.S2F4-mean(TT_sync_com.S2F4(zero_ref)),'Color',CusCol{2},'DisplayName','S_{ref,ch2}');
% Create lines for marking special areas of importance

% t_mark1=datetime("01.12.2022 06:07:06.000",'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone);
% t_mark2=datetime("02.12.2022 06:05:47.000",'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone);
% l1=line([t_mark1 t_mark1],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','DisplayName','');
% l2=line([t_mark2 t_mark2],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','DisplayName','');
% l1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% l2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Secondary axis for temperature

%plot(TT_sync_com.Time(idx.C25.dsc_1),TT_sync_com.S2F2(idx.C25.dsc_1)-mean(TT_sync_com.S2F2(zero_ref))+Offset_wl.S2F2,'DisplayName','C20-dsc_{cyc 1}','Color',cyan{4},'Marker','none','LineStyle','-','MarkerSize',3);
% plot(TT_sync_com.Time(idx.C25.chr_phase{3}),TT_sync_com.S2F2(idx.C25.chr_phase{3})-mean(TT_sync_com.S2F2(zero_ref))+Offset_wl.S2F2,'DisplayName','C20-chr4','Color',cyan{4},'Marker','none','LineStyle','-','MarkerSize',3);
%plot(TT_sync_com.Time(idx.C25.dsc),TT_sync_com.S2F2(idx.C25.dsc)-mean(TT_sync_com.S2F2(zero_ref))+Offset_wl.S2F2,'DisplayName','C20-dsc','Color',cyan{4},'Marker','none','LineStyle','-','MarkerSize',3);

% plot(TT_sync_com.Time(idx.C25.chr),TT_sync_com.S2F2(idx.C25.chr)-mean(TT_sync_com.S2F2(zero_ref))+Offset_wl.S2F2,'DisplayName','C20-DSC','Color',cyan{4},'Marker','none','LineStyle','-','MarkerSize',3);


yyaxis right
% plot(TT_sync_com.Time,TT_sync_com.("T_pico_mean_cell_01 / °C"),'DisplayName','T_{cell 1}','Color',[1 0.1 0.1],LineStyle='--');
plot(TT_sync_com.Time(1:867867),TT_sync_com.("T_pico_mean_cell_02 / °C")(1:867867),'DisplayName','T_{ambient hall}','Color',[1 0.1 0.1],LineStyle='--');  %GenI: Cell2 GenII : ambient hall
% plot(TT_sync_com.Time,TT_sync_com.("T1_pt100_degC"),'DisplayName','T_{pt100,T1}','Color',[1 0.1 0.8],LineStyle='-');
% plot(TT_sync_com.Time,TT_sync_com.T2,'DisplayName','T_{pt100,T2}','Color',[1 0.1 0.1],LineStyle='-');
%%%%%% commented June 16 by Luke %%%%%%%%%plot(TT_sync_com.Time,test,'DisplayName','T_{cali,T2}','Color',CusCol{2},LineStyle='--');
% plot(TT_sync_com.Time,TT_sync_com.("T_chamber_amb"),'DisplayName','T_{pt100,amb}','Color',[1 0.1 0.1],LineStyle='-'); %ambient chamber
plot(TT_sync_com.Time(1:867867),TT_sync_com.("T_pico_mean_ambient / °C")(1:867867),'DisplayName','T_{typeK,amb}','Color',[0.6 0.8 0.4],'LineStyle','--');


%     plot(synchronizedTimetable.Time,synchronizedTimetable.("Humidity"),'DisplayName','Humidity','Color',red{2},'LineWidth', 0.01); 
%     plot(synchronizedTimetable.Time,synchronizedTimetable.("y_cell_Temp"),'DisplayName','y_cell_Temp','Color',blue{2},'LineWidth', 0.01); 
%     plot(synchronizedTimetable.Time,synchronizedTimetable.("43_cell_Temp"),'DisplayName','43_cell_Temp','Color',blue{3},'LineWidth', 0.01); 
%     plot(synchronizedTimetable.Time,synchronizedTimetable.("x_cell_Temp"),'DisplayName','x_cell_Temp','Color',blue{4},'LineWidth', 0.01); 
%     plot(synchronizedTimetable.Time,synchronizedTimetable.("ambient_Temp"),'DisplayName','ambient_Temp','Color',margenta{5},'LineWidth', 0.01); 



ylabel('$\mathrm{T_{cell} \; / \; ^\circ{}C}$','Interpreter','latex')
legend('Location','east outside')
% xlim([TT_sync_com.Time(1) TT_sync_com.Time(end)]);

% ylim([24 34])
% xlim([TT_sync_com.Time(1) TT_sync_com.Time(end)]);



ax_handle.YAxis(2).Color=[1 0.1 0.1];

% Plot hero data in addtional subplot is available

if Hero_data
    ax2 = nexttile;
    ax_handle=gca;

    hold on;
    grid on;
    xlabel('$\mathrm{t_{exp} \; / \; HH:MM:SS}$','Interpreter','latex')%HH:MM:SS
    ylabel('$\mathrm{U \; / \; V}$','Interpreter','latex')
    plot(TT_sync_com.Time,TT_sync_com.("U / V(KO10)"),'DisplayName','U / V');

    % luke add:
    %scatter(TT_sync_com.Time( idx.C.chr_phase{1} ),TT_sync_com.("U / V")( idx.C.chr_phase{1} ),'MarkerEdgeColor', 'r'); %luke added 6 july

    yyaxis right

    plot(TT_sync_com.Time,TT_sync_com.("SOC / %(KO10)"),'DisplayName','SOC / %');
    legend('Location','eastoutside')
    
    linkaxes([ax1 ax2],'x');
    xlim([TT_sync_com.Time(1) TT_sync_com.Time(867867)])
%     xlim([t_low t_high])

end

%save
save_fig=0; %placeholder
if save_fig
    % exportname= input('What name should the saved figure have? Please write a string or char array as answer: ');
    ex_time=datetime(now,'ConvertFrom','datenum','TimeZone','Europe/Berlin','Format','ddMMyyyy_HHmmss');
    exportname_fbg=strcat(string(ex_time),"_FBG_Signals",".png");
    exportgraphics(fig_handle,fullfile(pic_folder_path,exportname_fbg),'Resolution',300)
end
% Alternative!
% Determine relevant vector beforehand and reduce lines drastically!

% % most general method:
% i=1;
% plot(TT_sync_com.Time, relevantDataTT.(datasetNames(i)) -mean(relevantDataTT.(datasetNames(i))(zero_ref)) ,'Color',margenta{4},'DisplayName',datasetNames(i),'LineWidth', 0.01,'Marker', 'none');
% 

% if Putty_data
%  
%     % plot selected favorites:
%     for i = 1:numel(datasetNames) 
%         if tempNames(i) == "DNE"
%             continue % if there is no temp data, skip
%         end
% 
%         plot(relevantDataTT.Time, relevantDataTT.(string(datasetNames(i) + "_T"  ) ) ,'DisplayName',string(datasetNames(i) + "_T"  ),'Color',colorList{i+1},'LineWidth', 0.01,'LineStyle', '-','Marker', '.'); 
%     end

%%additional code options from Luke
Luke_visualization_code % only inactive code for now
%% Analysis - Calibration
%% Strain calibration
%% 
% * Reduce number of function inputs?

if contains(Test_name,'C20') || contains(Test_name,'C25')
    strain_cal(Test_name,TT_sync_com,Meta.Kap,zero_ref,fig_siz,Offset_wl,margenta,cyan,save_fig,pic_folder_path,show_fit)
    % ^ Luke edited the end of the strain_cal script so that it works with his computer's path
end
% Identify Charge and Discharge Phases
% Added on 27 June 2023 by Luke McCarvill. Using findpeaks, identifies phases 
% and saves them in cell arrays ''|idx.C25.dsc_phase"| and "|idx.C25.chr_phase"| 
% based on the peak indices. (Or at least it did until I made it more general, 
% so it isn't just for the C25 data)

% run phaseIdentifier % right now this is just used to stop the script
%% Temperature calibration
% Indicate which sensors are used for temperature calibration (if applicable)
% Specify datasetNames and tempNames (temperature) 

if strcmp(Test_name,'20230906#25degC_CrateTest#ko10-xx_BLB1-xB') || strcmp(Test_name,'20230821#25degC_lowCRateTest#Ko10-xx_BLB1-xB')
    datasetNames = "S3F3"; % Edit as needed
    tempNames =    "x_cell_Temp"; % Edit as needed
elseif strcmp(Test_name,'20230906#25degC_CrateTest#ko10-xx_BLB1-xB')
    datasetNames = "S1F3"; % Edit as needed
    tempNames =    "x_cell_Temp"; % Edit as needed
else
    datasetNames = ["S1F2", "S2F1"        , "S3F3"       , "S4F4"       ]; % Edit as needed
    tempNames =    ["DNE" , "43_cell_Temp", "x_cell_Temp", "y_cell_Temp"]; % Edit as needed
    % IN THE FUTURE, MAKE THIS ^ AFFECT PuTTYsynchronizer
end
    % if there is no temp data available, mark it as "DNE" (does not
    % exist). These are the names from the PuTTYsynchronizer script.
if numel(datasetNames) ~= numel(tempNames)
    error('There must be a corresponding tempNames entry for each datasetNames entry. If a temperature vector is not available, indicate it as "DNE"')
end
%  - Now using Luke's temperature_cal_Putty()
% Using calibration technique is used here for the determination temperature 
% change 

if contains(Test_name,'Pyr') || contains(Test_name,'T_Ramp_BLB1_xB')  %make sure there is no pyr for electrical load
    temperature_cal_Putty(Test_name,TT_sync_com,zero_ref,colorList, colorListDark,save_fig,pic_folder_path,synchronizedTimetable,relevantDataTT,datasetNames,tempNames)
end

error("Stop here for now"); % stop the script here
% 
%% Visualization analysis results - work in progress
% Load calibration data
%% 
% * Make import from calibration data more clean
% * Save inverted calibration directly in calibration files - inversion of non 
% 1st order polynoms seems complicated
% * Right now: SOC from hero system - accurate enough calculated?
% * Right now SOC from charge is used during rest periods. Seems to be more 
% a hystersis depending on what was before.

if strcmp(Test_name,'20230102_25degC_1CC20_Zyk_Ko10-xx_S6_S3_S4')
%Temperature - T(wlc) - wlc=wavelengthchange
load('Y:\07_Experimente_und_Messdaten\Projekte\FiSens_PouchZelle\Kokam_10Ah\Analysis_results\Calibration\Temperature\20230109_25155015degC2Pyr-0C_Ko10-xx_S6_S3_S4_S1_T_cal.mat');

%invert calibration - wlc(T)
wlc_T=zeros(1,2);
wlc_T(1)=1/T_cal_pf.S2F2_1O(1);
wlc_T(2)=-T_cal_pf.S2F2_1O(2)/T_cal_pf.S2F2_1O(1);

%Strain
load('Y:\07_Experimente_und_Messdaten\Projekte\FiSens_PouchZelle\Kokam_10Ah\Analysis_results\Calibration\Strain\20230102_25degC_1CC20_Zyk_Ko10-xx_S6_S3_S4_strain_cal.mat');


%%Use calibration to determine signal(wlc) from temperature and state of charge

% T influce
T_cal_ref=mean(TT_sync_com.T2(zero_ref));
T_change_vek=TT_sync_com.T2;
wlc_T_change=polyval(wlc_T,T_change_vek);

% SOC influence
wlc_SOC_change=zeros(length(TT_sync_com.("SOC / %")),1);
idx_comp.chr=find(TT_sync_com.("I / A")<=0);
idx_comp.dsc=find(TT_sync_com.("I / A")>0);
wlc_SOC_change(idx_comp.chr)=polyval(wlc_strain_cal_pf.S2F2_chr_1O,TT_sync_com.("SOC / %")(idx_comp.chr)/100);
wlc_SOC_change(idx_comp.dsc)=polyval(wlc_strain_cal_pf.S2F2_dsc_7O,TT_sync_com.("SOC / %")(idx_comp.dsc)/100);

% Determine temperature after complete calibration
test=polyval(T_cal_pf.S2F2_1O,TT_sync_com.S2F2-mean(TT_sync_com.S2F2(zero_ref))-wlc_SOC_change);
end
% Use calibration data for signal contribution analysis
% With calibrations from the same setup, the signal contributions should be 
% completly describable.
%% 
% * Does it work if calibration and measurement are not done in same consecutive 
% measurement - maybe some kind of hysteresis?
% Create temperature curve over experiment
% Determine also the fit of the calibrated measurements to the reference measurements

[T_SP1_cu_1thO_corr,sig_SP1_cu_SP1O_corr]=polyval(T_SP1_fit_para_1thO_corr,TT_sync_com.S3F1-mean(TT_sync_com.S3F1(zero_ref))-0.25*Ambient_fixed_corr.S3,err_SP1_1O_corr);


figure()
fig_handle=gcf;
set(gcf,'Visible','on')
set(gcf,'position',fig_siz)
set(gcf,'color','w');
set(0, 'DefaultLineLineWidth', 2.5);
set(gcf,'DefaultAxesFontSize',GCAfontsize);
hold on;
grid on;
ax_handle=gca;


xlabel('$\mathrm{t_{exp} \; / \; HH:MM:SS}$','Interpreter','latex')
ylabel('$\mathrm{T \; / \; ^{\circ}C}$','Interpreter','latex')
plot(TT_sync_com.Time,TT_sync_com.T3,'r','DisplayName','pt100-ref');
plot(TT_sync_com.Time,T_SP1_cu_1thO_corr,'k','DisplayName','S4P-lin-cr-fit');
% shadedErrorBar(TT_sync_com.Time,T_SP1_cu_1thO_corr,sig_SP1_cu_SP1O_corr,'lineProps', {'DisplayName','LinFit-S4P_{T,amb}','Color',red{5},'LineWidth',2.5})%reference TO and lambda0 of calibration // +T_ref_pt100
% shadedErrorBar läuft so nicht mit der Datenmenge anscheind
xlim([TT_sync_com.Time(1) TT_sync_com.Time(end)])
legend('Location','best')

figure()
fig_handle=gcf;
set(gcf,'Visible','on')
set(gcf,'position',fig_siz)
set(gcf,'color','w');
set(0, 'DefaultLineLineWidth', 2.5);
set(gcf,'DefaultAxesFontSize',GCAfontsize);
hold on;
grid on;
ax_handle=gca;

xlabel('$\mathrm{t_{exp} \; / \; HH:MM:SS}$','Interpreter','latex')
ylabel('$\Delta\mathrm{T \; / \; ^{\circ}C}$','Interpreter','latex')
plot(TT_sync_com.Time,TT_sync_com.T3-T_SP1_cu_1thO_corr,'DisplayName','T_{ref}-T_{S4P,cr}','Color',green{4});
plot(TT_sync_com.Time,TT_sync_com.T3-polyval(T_SP1_fit_para_1thO_corr,TT_sync_com.S3F1-mean(TT_sync_com.S3F1(zero_ref)),err_SP1_1O_corr),'DisplayName','T_{ref}-T_{S4P}','Color',cyan{4});
legend('Location','best')
ylim([-0.8 0.8])
xlim([TT_sync_com.Time(1) TT_sync_com.Time(end)])

yyaxis right
ylabel('T_{hall,diff}')
ylim([-5 5])
plot(TT_sync_com.Time,TT_sync_com.("T_pico_mean_cell_02 / °C")-mean(TT_sync_com.("T_pico_mean_cell_02 / °C")(zero_ref)),'DisplayName','T_{hall,amb}','Color',[1 0.1 0.1],LineStyle='--');  % ambient hall
ax_handle.YAxis(2).Color=[1 0.1 0.1];
% Claibration test for specific temperature

bla % only for simple stop of code
if strcmp(Test_name,'15dyn50degC-0C_Ko10-xx_F1S5GenII_20221212')
    load('.\ParameterSets\Temp_calibration')
    hall_phase.idx=1801630:length(TT_sync_com.Time);
    hall_phase.time=TT_sync_com.Time(hall_phase.idx);
    hall_phase.Tamb_hall=TT_sync_com.("T_pico_mean_ambient / °C")(hall_phase.idx);
    hall_phase.Tamb_chamber=TT_sync_com.T3(hall_phase.idx);
    hall_phase.SP4=TT_sync_com.S3F1(hall_phase.idx)-mean(TT_sync_com.S3F1(zero_ref));

    T0_hall=18.5; % average of the daytemperature of test-start
    hall_inf=1.8e-3; % pm/K
    hall_phase.hall_change = -(hall_phase.Tamb_hall-T0_hall)*hall_inf+hall_phase.SP4;

    %% fit experiment for dynamic phase
    t_fit_low=("08.12.2022 07:24:57.801");
    t_fit_high=("10.12.2022 08:58:32.777");
    idx_fit_1stpeak=find(TT_sync_com.Time==datetime(t_fit_low,'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone));
    idx_fit_end=find(TT_sync_com.Time==datetime(t_fit_high,'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone));

    fit_lambda_x_test=TT_sync_com.S3F1(idx_fit_1stpeak:idx_fit_end)-mean(TT_sync_com.S3F1(zero_ref));
    fit_T_y_test=TT_sync_com.T3(idx_fit_1stpeak:idx_fit_end)-mean(TT_sync_com.T3(zero_ref));

    T_cal_fit_6thO_test=polyfit(TT_sync_com.S3F1(idx_fit_1stpeak:idx_fit_end)-mean(TT_sync_com.S3F1(zero_ref)),TT_sync_com.T3(idx_fit_1stpeak:idx_fit_end)-mean(TT_sync_com.T3(zero_ref)),5);
    %%

    figure()
    fig_handle=gcf;
    set(gcf,'Visible','on')
    set(gcf,'position',fig_siz)
    set(gcf,'color','w');
    set(0, 'DefaultLineLineWidth', 2.5);
    set(gcf,'DefaultAxesFontSize',GCAfontsize);
    hold on;
    grid on;
    plot(TT_sync_com.S3F2(idx_fit_1stpeak:idx_fit_end)-mean(TT_sync_com.S3F2(zero_ref)),TT_sync_com.T3(idx_fit_1stpeak:idx_fit_end),'DisplayName','dyn deviation');
%     plot(hall_phase.SP4,hall_phase.Tamb_chamber,'DisplayName','corrected deviation');
%     plot(hall_phase.hall_change,hall_phase.Tamb_chamber,'DisplayName','corrected deviation');
    
    figure()
    fig_handle=gcf;
    set(gcf,'Visible','on')
    set(gcf,'position',fig_siz)
    set(gcf,'color','w');
    set(0, 'DefaultLineLineWidth', 2.5);
    set(gcf,'DefaultAxesFontSize',GCAfontsize);
    hold on;
    ax_handle=gca;
    grid on;
    xlabel('$\mathrm{t_{exp} \; / \; HH:MM:SS}$','Interpreter','latex')
    ylabel('$\Delta\mathrm{\lambda \; / \; nm}$','Interpreter','latex')
    plot(hall_phase.time,hall_phase.SP4,'Color',[1 0.6 0.1],'DisplayName','T_{S4P,amb}');
    plot(hall_phase.time,hall_phase.hall_change,'Color',[1 0.6 0.9],'DisplayName','T_{S4P,amb,hc}');

    yyaxis right
    plot(hall_phase.time,hall_phase.Tamb_chamber,'DisplayName','T_{pt100,amb}','Color',[1 0.1 0.1]);
    plot(hall_phase.time,hall_phase.Tamb_hall,'DisplayName','T_{typeK,amb}','Color',[1 0.4 0.4],'LineStyle','--');
    ylabel('$\mathrm{T_{cell} \; / \; ^\circ{}C}$','Interpreter','latex')
    legend('Location','best')

    figure()
    fig_handle=gcf;
    set(gcf,'Visible','on')
    set(gcf,'position',fig_siz)
    set(gcf,'color','w');
    set(0, 'DefaultLineLineWidth', 2.5);
    set(gcf,'DefaultAxesFontSize',GCAfontsize);
    hold on;
    ax_handle=gca;
    plot(TT_sync_com.Time,TT_sync_com.T3);
    plot(TT_sync_com.Time,polyval(T_SP4_fit_para_6thO,TT_sync_com.S3F1-lambda_ref_S3F1)+T_ref_pt100); %reference TO and lambda0 of calibration
    plot(TT_sync_com.Time,polyval(T_SP4_fit_para_6thO,TT_sync_com.S3F1-mean(TT_sync_com.S3F1(zero_ref)))+mean(TT_sync_com.T3(zero_ref))); %reference of measurement
    plot(TT_sync_com.Time,polyval(T_cal_fit_6thO_test,TT_sync_com.S3F1-mean(TT_sync_com.S3F1(zero_ref)))+mean(TT_sync_com.T3(zero_ref))); %reference of measurement and calibration for test
end
bla
% 
% Create relative time vector for better visualization on poster or presentation
% TODO:
%% 
% * Rework und don't hardcode times

t_lim_low=datetime("30.10.2022 18:45:00.000",'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone);
t_lim_mid_alt=datetime("26.09.2022 08:43:11.000",'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone);
t_lim_mid=datetime("26.09.2022 08:44:13.000",'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone);
t_lim_high=datetime("01.11.2022 10:42:00.000",'InputFormat','dd.MM.yyyy HH:mm:ss.SSS','Format',['dd.MM.yyyy HH:mm:ss.SSS'],'TimeZone',TT_sync_com.Time.TimeZone);
% xlim([t_lim_low,t_lim_high])

%save

% exportname= input('What name should the saved figure have? Please write a string or char array as answer: ');
% exportgraphics(fig_handle,fullfile('Z:\06 Modellierung\Matlab\Analysis_scripts\Analysis_pics\PouchCell',exportname,'.png'),'Resolution',300)
%% 
% Show only part of complete profile

Cyc_part.idx=find(TT_sync_com.Time>=t_lim_low & TT_sync_com.Time<=t_lim_high);
Cyc_part.numel=numel(Cyc_part.idx);
Cyc_part.t_max=(datenum(t_lim_high)-datenum(t_lim_low))*24;
Cyc_part.t_dur=0:Cyc_part.t_max/(Cyc_part.numel-1):Cyc_part.t_max;

fig=figure("Name",'FBG wavelength shift part','DefaultAxesFontSize',GCAfontsize,'DefaultAxesFontName','Verdana');
fig_handle=gcf;


% gca.fontsize=GCAfontsize;
%Voltage and current during 0.5 C charging phase and break
tiledlayout(2,1)
set(gcf,'Visible','on')
set(gcf,'position',fig_siz)
set(gcf,'color','w');
set(0, 'DefaultLineLineWidth', 2);
%set(gcf,'DefaultAxesFontSize',22);

ax1 = nexttile;
ax_handle=gca;
hold on;
grid on;
xlabel('$\mathrm{t \; / \; h}$','Interpreter','latex')
ylabel('$\Delta\mathrm{\lambda \; / \; nm}$','Interpreter','latex')

%plot(TT_sync_com.Time,TT_sync_com.S1F1-mean(TT_sync_com.S1F1(zero_ref)));
% plot(TT_sync_com.Time,TT_sync_com.S1F2-mean(TT_sync_com.S1F2(zero_ref)));
% plot(TT_sync_com.Time,TT_sync_com.S1F3-mean(TT_sync_com.S1F3(zero_ref)),'Color','b');
% plot(TT_sync_com.Time,TT_sync_com.S1F4-mean(TT_sync_com.S1F4(zero_ref)));
% plot(TT_sync_com.Time,TT_sync_com.S1F5-mean(TT_sync_com.S1F5(zero_ref)),'Color',CusCol{3});
plot(Cyc_part.t_dur,TT_sync_com.S1F6(Cyc_part.idx)-mean(TT_sync_com.S1F6(zero_ref)),'DisplayName','$\mathrm{S1_{T, \varepsilon}}$','Color',CusCol{1});
plot(Cyc_part.t_dur,TT_sync_com.S1F5(Cyc_part.idx)-mean(TT_sync_com.S1F5(zero_ref)),'DisplayName','$\mathrm{S2_{T}}$','Color',CusCol{2});
plot(Cyc_part.t_dur,ficorr.S1F6(Cyc_part.idx)-mean(TT_sync_com.S1F6(zero_ref)),'DisplayName','$\mathrm{S1_{\epsilon}}$','Color',CusCol{4});
% plot(TT_sync_com.Time,TT_sync_com.S2F1-mean(TT_sync_com.S2F1(zero_ref)),'Color',Blue{1});
% plot(TT_sync_com.Time,TT_sync_com.S1F2-mean(TT_sync_com.S1F2(zero_ref)));
% plot(TT_sync_com.Time,TT_sync_com.W9-mean(TT_sync_com.W9(zero_ref)),'Color',Blue{2});
% plot(TT_sync_com.Time,TT_sync_com.S1F4-mean(TT_sync_com.S1F4(zero_ref)));
% plot(TT_sync_com.Time,TT_sync_com.S2F5-mean(TT_sync_com.S2F5(zero_ref)),'Color',Blue{3});
% plot(TT_sync_com.Time,TT_sync_com.S2F6-mean(TT_sync_com.S2F6(zero_ref)),'Color',Blue{3});
ylim([-0.05 0.15])

yyaxis right
plot(Cyc_part.t_dur,TT_sync_com.("T_pico_mean_cell_01 / °C")(Cyc_part.idx),'DisplayName','$\mathrm{T_{cell}}$','Color',CusCol{3});
% plot(Cyc_part.t_dur,Pico_TT.("T_pico_mean_ambient / °C")(Cyc_part.idx),'DisplayName','T_{amb}');
ylabel('$\mathrm{T_{cell} \; / \; ^\circ{}C}$','Interpreter','latex')
legend('Location','best','Interpreter','latex')
ylim([25 31])
%legend('FBG_{mid} - T_{ref}','FBG_{mid} - \epsilon, T', 'Location','eastoutside','NumColumns',1) %,'FBG_{mid} - \epsilon'
%legend('T_{top,right}','T_{middle}','T_{left}','T_{bot}')
ax_handle.YAxis(2).Color=CusCol{3};

% Create background colors
% The f inding algorithm needs some work. Use find?

% eps=0.5;
% [JP]=ischange(TT_sync_com.("I / A")(Cyc_part.idx),'linear','Threshold',5);
% POI=Cyc_part.t_dur(JP);
% Disc_CC_start=POI(2);
% Disc_CC_end=POI(3);
% Char_CC_start=POI(3);
% Char_CC_end=Cyc_part.t_dur(end);
% text((Disc_CC_start+Disc_CC_end)/2,0.99*ax_handle.YLim(2),'Dsc CC','HorizontalAlignment','center','FontSize',GCAfontsize);
% text((Char_CC_start+Char_CC_end)/2,0.99*ax_handle.YLim(2),'Chr CC','HorizontalAlignment','center','FontSize',GCAfontsize);
% fill([Disc_CC_start,Disc_CC_start,Disc_CC_end,Disc_CC_end],[ax_handle.YLim(2),ax_handle.YLim(1),ax_handle.YLim(1),ax_handle.YLim(2)],'r','Edgealpha',0.08,'Facealpha',0.08,'HandleVisibility','off');
% fill([Char_CC_start,Char_CC_start,Char_CC_end,Char_CC_end],[ax_handle.YLim(2),ax_handle.YLim(1),ax_handle.YLim(1),ax_handle.YLim(2)],'g','Edgealpha',0.08,'Facealpha',0.08,'HandleVisibility','off');

% Curves of electrical quanitites

ax2 = nexttile;
ax_handle2=gca;
hold on;
grid on;
xlabel('$\mathrm{t \; / \; h}$','Interpreter','latex')
ylabel('$\mathrm{I \; / \; A}$','Interpreter','latex')
plot(Cyc_part.t_dur,TT_sync_com.("I / A")(Cyc_part.idx),'DisplayName','I / A','Color',TuCol{1});
yyaxis right
plot(Cyc_part.t_dur,TT_sync_com.("U / V")(Cyc_part.idx),'DisplayName','U / V','Color',TuCol{2});
ylabel('$\mathrm{U \; / \; V}$','Interpreter','latex')
legend('Location','eastoutside')
ax_handle2.YAxis(2).Color=TuCol{2}
linkaxes([ax1 ax2],'x')
xlim([0 Cyc_part.t_max])
% ylim([-10.1 10.1])

% 
% Temperatureprofile
% Sort by ambient temperatur for having same conditions. Investigate the effect 
% drift during different phases during the pyramide steps.

% set parameters
T_amb_set=[15 20 25 30 35 40];
eps_T=0.5; % increment of viable temperature range for comparison

T_amb_val=cell(1,length(T_amb_set));
T_amb_idx=cell(1,length(T_amb_set));
T_amb_varnames=cell(1,length(T_amb_set));

for i=1:length(T_amb_set)
    T_amb_varnames{1,i}=['deg_',num2str(T_amb_set(i)),'_C'];
    sprintf('Temperature %i °C',T_amb_set(i));
    T_amb_find=find(TT_sync_com.("T_pico_mean_ambient / °C") <= T_amb_set(i)+eps_T & TT_sync_com.("T_pico_mean_ambient / °C") >= T_amb_set(i)-eps_T);
    T_amb_idx{:,i}=T_amb_find;
    T_amb_val{:,i}=TT_sync_com.("T_pico_mean_ambient / °C")(T_amb_find);
%put into table
    if i==length(T_amb_set)
        T_amb_idx_table=cell2table(T_amb_idx(:,:),'VariableNames',T_amb_varnames);
        T_amb_val_table=cell2table(T_amb_val(:,:),'VariableNames',T_amb_varnames);
    end
end
%% 
% Which temperature would you like to look at in detail?
% 
% Time shif at 25 °C because one logging started earlier then the other...catch 
% error

k=3;
clear switch_pl switch_idx P1 P2 P3 P
%t_rel_T=0:1:numel(T_amb_idx_table.(k){1,1})-1;
[switch_pl,switch_idx]=findpeaks(diff(T_amb_idx_table.(k){1,1}),'MinPeakHeight',10e3);


% Determine periods of time of interest
Before_switch=T_amb_idx_table.(k){1,1}(switch_idx);%T_amb_idx_table.(k){1,1}(1)+switch_idx
After_switch=T_amb_idx_table.(k){1,1}(switch_idx+1);%T_amb_idx_table.(k){1,1}(1)+switch_idx;

if numel(switch_idx)==2
    P1.idx=T_amb_idx_table.(k){1,1}(1):Before_switch(1);
    P2.idx=After_switch(1):Before_switch(2);
    P3.idx=After_switch(2):T_amb_idx_table.(k){1,1}(end);
    P.idx=[P1.idx P2.idx P3.idx];
    P1.Time=TT_sync_com.Time(P1.idx); %cooling
    P2.Time=TT_sync_com.Time(P2.idx); %heating
    P3.Time=TT_sync_com.Time(P3.idx); %cooling
    P1.trel=[0,cumsum(diff(datenum(P1.Time)))'];
    P2.trel=[0,cumsum(diff(datenum(P2.Time)))'];
    P3.trel=[0,cumsum(diff(datenum(P3.Time)))'];
    P.trel=[P1.trel';P2.trel'+P1.trel(end);P3.trel'+P2.trel(end)+P1.trel(end)]*24; %*60*60
else
    P1.idx=T_amb_idx_table.(k){1,1}(1):Before_switch(1); %heating
    P2.idx=After_switch(1):T_amb_idx_table.(k){1,1}(end); %cooling
    P.idx=[P1.idx P2.idx];
    P1.Time=TT_sync_com.Time(P1.idx); %cooling
    P2.Time=TT_sync_com.Time(P2.idx); %heating
    P1.trel=[0,cumsum(diff(datenum(P1.Time)))'];
    P2.trel=[0,cumsum(diff(datenum(P2.Time)))'];
    P.trel=[P1.trel';P2.trel'+P1.trel(end)]*24; %*60*60
end
% Show temperature profiles and signals from fibre optical sensors

fig=figure("Name",'FBG wavelength shift','DefaultAxesFontSize',GCAfontsize);
set(gcf,'color','w');
fig_handle=gcf;
set(0, 'DefaultLineLineWidth', 1.5);
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'position',[0 0 2560 1440])
% gca.fontsize=GCAfontsize;
%Voltage and current during 0.5 C charging phase and break
tiledlayout(3,1)
set(gcf,'Visible','on')


ax1 = nexttile;
ax_handle=gca;
hold on;
grid on;
%xlabel('$\mathrm{t_{exp} \; / \; HH:MM:SS}$','Interpreter','latex')
% ylabel('$\Delta\mathrm{\lambda \; / \; nm}$','Interpreter','latex')

%xlim([duration([0,0,0]) duration([110,0,0])])
%plot(FiSens_spectra_Ana.Time,FiSens_spectra_Ana{:,1:11}-FiSens_spectra_Ana{1,1:11})

plot(P.trel,TT_sync_com.("T_pico_mean_ambient / °C")(P.idx));
plot(P.trel,TT_sync_com.('T_pico_mean_cell_01 / °C')(P.idx))
%plot(t_rel_T/3600,Pico_TT.(2)(T_amb_idx_table.(k){1,1}))
xlabel('$\mathrm{t \; / \; h}$','Interpreter','latex')
ylabel('$\mathrm{T \; / \; ^\circ{}C}$','Interpreter','latex')
%ylabel('$\mathrm{T_{cell} \; / \; ^\circ{}C}$','Interpreter','latex')
% xlim([t_rel_T(1)/3600 t_rel_T(end)/3600])
if length(switch_idx)==2
    line([P.trel(switch_idx(1)) P.trel(switch_idx(1))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--')
    line([P.trel(switch_idx(2)) P.trel(switch_idx(2))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--')
    legend('amb','cell_1','','','','Location','best')
else
    line([P.trel(switch_idx(1)) P.trel(switch_idx(1))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--')
    legend('amb','cell_1','cell_2','','Location','best')
end

ax2 = nexttile;
ax_handle=gca;
hold on;
grid on;
xlabel('$\mathrm{t \; / \; h}$','Interpreter','latex')
ylabel('$\Delta\mathrm{\lambda \; / \; nm}$','Interpreter','latex')
plot(P.trel,TT_sync_com.S1F6(P.idx)-mean(TT_sync_com.S1F6(zero_ref)),'DisplayName','S1_{T, \epsilon}','Color','b','LineStyle','--');
% plot(P.trel,TT_sync_com.S1F5(P.idx)-mean(TT_sync_com.S1F5(zero_ref)),'DisplayName','S2_{T}','Color','k');
plot(P.trel,ficorr.S1F6(P.idx)-mean(TT_sync_com.S1F6(zero_ref)),'DisplayName','S1_\epsilon','Color','c');
% plot(P.trel,TT_sync_com.S1F4(P.idx)-mean(TT_sync_com.S1F4(zero_ref)),'DisplayName','S3');
% plot(P.trel,TT_sync_com.S1F3(P.idx)-mean(TT_sync_com.S1F3(zero_ref)),'DisplayName','S4');
% plot(P.trel,ficorr.S1F2(P.idx)-mean(TT_sync_com.S1F2(zero_ref)),'DisplayName','S5_\epsilon','Color','k');
% plot(P.trel,TT_sync_com.S1F2(P.idx)-mean(TT_sync_com.S1F2(zero_ref)),'DisplayName','S5','Color','k','LineStyle','--');
% plot(P.trel,TT_sync_com.S1F1(P.idx)-mean(TT_sync_com.S1F1(zero_ref)),'DisplayName','S6');
if length(switch_idx)==2
    legend('Location','eastoutside')
    line([P.trel(switch_idx(1)) P.trel(switch_idx(1))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','HandleVisibility','off')
    line([P.trel(switch_idx(2)) P.trel(switch_idx(2))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','HandleVisibility','off')
else
    legend('Location','eastoutside')
    line([P.trel(switch_idx(1)) P.trel(switch_idx(1))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','HandleVisibility','off')
end

% xlim([t_rel(1) t_rel(9000000)])

ax3 = nexttile;
hold on;
grid on;
xlabel('$\mathrm{t \; / \; h}$','Interpreter','latex')
ylabel('$\mathrm{I \; / \; A}$','Interpreter','latex')
plot(P.trel,TT_sync_com.("I / A")(P.idx));
if length(switch_idx)==2
    line([P.trel(switch_idx(1)) P.trel(switch_idx(1))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','DisplayName','')
    line([P.trel(switch_idx(2)) P.trel(switch_idx(2))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','DisplayName','')
else
    line([P.trel(switch_idx(1)) P.trel(switch_idx(1))],[ax_handle.YLim(1) ax_handle.YLim(2)],'Color','black','LineStyle','--','DisplayName','')
end

yyaxis right
plot(P.trel,TT_sync_com.("U / V")(P.idx));
ylabel('$\mathrm{U \; / \; V}$','Interpreter','latex')

linkaxes([ax1 ax2 ax3],'x')
xlim([P.trel(1) P.trel(end)])
% Show other temperature sensor

% if numfib

%