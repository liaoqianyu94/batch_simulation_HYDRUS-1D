clc
clear all

work_dir='F:\HYDRUS-1D\work_dir\';
file_name='run_1';
file_path=strcat(work_dir,file_name)
exe_path='F:\HYDRUS-1D\H1D_CALC.EXE';
mkdir(file_path);

V_inform_path='F:\simulation_ameriflux\V_data.xlsx';
site='US-Ne3';
year=2007;
[p_height,p_LAI,p_RD,year_span]=read_V(V_inform_path,site,year);

obs_csv_path='F:\simulation_ameriflux\AMF_US-Ne3_BASE-BADM_13-5\AMF_US-Ne3_BASE_HR_13-5.csv';
simulation_strart_time=year_span(1)*10000;
simulation_end_time=20071101*10000;
[BC,MC,albedo,sm,ET,ET_balance]=read_csv(obs_csv_path,simulation_strart_time,simulation_end_time);
% p_RD=-100*(max(p_RD)-p_RD)./(max(p_RD)-min(p_RD));
p_RD=-1*p_RD-4;
p_RD(find(p_RD<0))=0;
albedo(find(isnan(albedo)))=0.23;

% nsm=0.1*sm(:,1)+0.15*sm(:,2)+0.25*sm(:,3)+0.5*sm(:,4);
% for i=1:length(nsm)-1
%     diff_sm(i)=nsm(i+1)-nsm(i)-0.01;
% end
%     diff_sm(find(diff_sm < 0))=0;
%     diff_sm=[diff_sm,0];
%     diff_sm(find(isnan(diff_sm)))=0;
%     BC(:,2)=diff_sm;

MC=[MC,p_height,albedo,p_LAI,p_RD*(-1)];
hour=num2str(length(p_RD));

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% write selector.in
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 file_ID=fopen(strcat(file_path,'\','SELECTOR.IN'),'wt');
fprintf(file_ID,'Pcp_File_Version=4\n');
fprintf(file_ID,'*** BLOCK A: BASIC INFORMATION *****************************************\n');
fprintf(file_ID,'Heading\n');
fprintf(file_ID,'Simulating root soil moisture\n');
fprintf(file_ID,'LUnit  TUnit  MUnit  (indicated units are obligatory for all input data)\n');
fprintf(file_ID,'cm\n'); % Length unit
fprintf(file_ID,'hours\n'); % Time unit
fprintf(file_ID,'mmol\n'); % Mass unit for concentration
fprintf(file_ID,'lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse\n');
fprintf(file_ID,' t     f     f      t     t     f      f     f       t         t         f\n');
fprintf(file_ID,'lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy\n');
fprintf(file_ID,' f       f       t       t       f       f       f       f       f       f\n');
fprintf(file_ID,'NMat    NLay  CosAlpha\n');
fprintf(file_ID,'  1       1       1\n');

fprintf(file_ID,'*** BLOCK B: WATER FLOW INFORMATION ************************************\n');
fprintf(file_ID,'MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)\n');
fprintf(file_ID,'  50    0.001      1\n');
fprintf(file_ID,'TopInf WLayer KodTop InitCond\n');
fprintf(file_ID,' t     t      -1       t\n');
fprintf(file_ID,'BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep\n');
fprintf(file_ID,' f     f     t     f     -1      f      0\n');
fprintf(file_ID,'    hTab1   hTabN\n');
fprintf(file_ID,'    1e-006   10000\n');
fprintf(file_ID,'    Model   Hysteresis\n');
fprintf(file_ID,'      0          0\n');
fprintf(file_ID,'   thr     ths    Alfa      n         Ks       l\n');
fprintf(file_ID,'  0.143    0.49   0.05    1.51     2.16   6.98\n');
fprintf(file_ID,'*** BLOCK C: TIME INFORMATION ******************************************\n');
fprintf(file_ID,'        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL\n');
fprintf(file_ID,'        0.01      0.003        1     1.3     0.7     3     7     1\n');
fprintf(file_ID,'      tInit        tMax\n');
fprintf(file_ID,['          0        ',hour,'\n']);
fprintf(file_ID,'  lPrintD  nPrintSteps tPrintInterval lEnter\n');
fprintf(file_ID,'     f           1            1       f\n');
fprintf(file_ID,'TPrint(1),TPrint(2),...,TPrint(MPL)\n');
fprintf(file_ID,['       ',hour,'\n']); % TPrint: First specified print-time

fprintf(file_ID,'*** BLOCK D: ROOT GROWTH INFORMATION ***********************************\n');
fprintf(file_ID,'iRootDepthEntry\n'); 
fprintf(file_ID,'        2\n');
fprintf(file_ID,'     iRFak     tRMin     tRMed     tRMax     xRMin     xRMed     xRMax   tPeriod\n');
fprintf(file_ID,'        0      265       3000     4392      1        120.9       121     4127\n');

fprintf(file_ID,'*** BLOCK G: ROOT WATER UPTAKE INFORMATION *****************************\n');
fprintf(file_ID,'     Model  (0 - Feddes, 1 - S shape)  cRootMax    OmegaC\n');
fprintf(file_ID,'        0                                   1\n');
fprintf(file_ID,'       P0       P2H       P2L       P3          r2H        r2L\n');
% fprintf(file_ID,'      -10      -750     -2000    -16000   0.0208333  0.00416667\n');
% fprintf(file_ID,'POptm(1),POptm(2),...,POptm(NMat)\n');
% fprintf(file_ID,'     -25\n');
fprintf(file_ID,'      -15      -500      -500     -8000   0.5  0.1\n');
fprintf(file_ID,'POptm(1),POptm(2),...,POptm(NMat)\n');
fprintf(file_ID,'     -30\n');

fprintf(file_ID,"*** END OF INPUT FILE 'SELECTOR.IN' ************************************");
fclose(file_ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% write atmosph.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_ID=fopen(strcat(file_path,'\','ATMOSPH.IN'),'wt'); 
fprintf(file_ID,'Pcp_File_Version=4\n');
fprintf(file_ID,'*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************\n');
fprintf(file_ID,'   MaxAL                    (MaxAL = number of atmospheric data-records)\n');
fprintf(file_ID,['   ',hour,'\n']);
fprintf(file_ID,'DailyVar  SinusVar  lLay  lBCCycles lInterc lDummy  lDummy  lDummy  lDummy  lDummy\n');
fprintf(file_ID,'       f       f       f       f       t       f       f       f       f       f\n');
fprintf(file_ID,' aIntercep (0.25 mm)\n');
fprintf(file_ID,'    0.25\n');
fprintf(file_ID,' hCritS                 (max. allowed pressure head at the soil surface)\n');
fprintf(file_ID,'      0\n');
fprintf(file_ID,'       tAtm        Prec       rSoil       rRoot      hCritA          rB          hB          ht        tTop        tBot        Ampl    RootDepth\n');
for i=1:size(BC,1)
fprintf(file_ID,'%12d %12.3f %12d %12d %12d %12d %12d %12d\n',BC(i,:));
end
fprintf(file_ID,"end*** END OF INPUT FILE 'ATMOSPH.IN' **********************************");
fclose(file_ID);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% write meteo.in
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 file_ID=fopen(strcat(file_path,'\','METEO.IN'),'wt');
fprintf(file_ID,'Pcp_File_Version=4\n');
fprintf(file_ID,'* METEOROLOGICAL PARAMETERS AND INFORMATION |||||||||||||||||||||||||||||||\n');
fprintf(file_ID,' MeteoRecords Radiation Penman-Hargreaves\n');
fprintf(file_ID,['         ',hour,'        2       f\n']);
fprintf(file_ID,'  lEnBal  lDaily  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy\n');
fprintf(file_ID,'       f       f       f       f       f       t       f       f       f       f\n');
fprintf(file_ID,' WindHeight     TempHeight\n');
fprintf(file_ID,'        550            550\n');
fprintf(file_ID,' iCrop (=0: no crop, =1: constant, =2: table, =3: daily)  SunShine  RelativeHum\n');
fprintf(file_ID,'         3                                                3         0\n');
fprintf(file_ID,'iLai (=0: given, =1: grass, =2; alfalfa, =3: surface fraction)  rExtinct\n');
fprintf(file_ID,'    0                                                              0.463\n');
fprintf(file_ID,' Interception\n');
fprintf(file_ID,'    1\n');
fprintf(file_ID,'aInterc [mm]\n');
fprintf(file_ID,'      0.25\n');
fprintf(file_ID,'Daily values\n');
fprintf(file_ID,'       t        Rad        TMax        TMin     RHMean      Wind    SunHours CropHeight     Albedo   LAI(SCF)      rRoot\n');
fprintf(file_ID,'      [T]  [MJ/m2/d]       [C]         [C]       [%%]     [km/d]     [hour]      [L]           [-]        [-]        [L]\n');
for i=1:size(MC,1)
fprintf(file_ID,'%12d %12.3f %12.2f %12.2f %12.1f %12.1f %12d %12.2f %12.2f %12.1f %12.1f\n',MC(i,:));
end
fprintf(file_ID,"end *** END OF INPUT FILE 'METEO.IN' **********************************\n");
fclose(file_ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% write profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%

nodel_n=101;
nodel_number=1:1:nodel_n;nodel_number=nodel_number';
x_n=0:-1:-nodel_n+1;x_n=x_n';
first_sm=[repelem(sm(1,1)./100,10,1);repelem(sm(1,2)./100,15,1);...
         repelem(sm(1,3)./100,25,1);repelem(sm(1,4)./100,51,1)];

Initial_sm=first_sm;
Mat=repelem(1,nodel_n,1);
Lay=repelem(1,nodel_n,1);
Beta=[repelem(0,5,1);repelem(1,nodel_n-5,1)];
% x=linspace(200,0,201);x=x';
% L=linspace(0,200,201);L=L';
% Lr=125;
% 
% for i=1:201
% if x(i)>L(i)-0.2*Lr
%     Beta(i,1)=1.667/Lr;
% elseif x(i)<=L(i)-0.2*Lr &  x(i)>=L(i)-Lr
%     Beta(i,1)=(2.0833/Lr)*(1-(L(i)-x(i))./Lr);
% elseif x(i)<L(i)-Lr
%     Beta(i,1)=0;
% end
% end

Axz=repelem(1,nodel_n,1);
Bxz=repelem(1,nodel_n,1);
Dxz=repelem(1,nodel_n,1);
Tem=repelem(10,nodel_n,1);

SPS=[nodel_number,x_n,Initial_sm,Mat,Lay,Beta,Axz,Bxz,Dxz];


file_ID=fopen(strcat(file_path,'\','PROFILE.DAT'),'wt');
fprintf(file_ID,'Pcp_File_Version=4\n');
fprintf(file_ID,'    2\n');% Number of fixed nodes.
fprintf(file_ID,'    1  0.000000e+000  1.000000e+000  1.000000e+000\n'); 
fprintf(file_ID,'    2 -1.000000e+002  1.000000e+000  1.000000e+000\n');
fprintf(file_ID,'  101    0    0    1 x         h      Mat  Lay      Beta           Axz            Bxz            Dxz          Temp          Conc \n');
for i=1:size(SPS,1)
fprintf(file_ID,'%8d %8d %8.2f %8d %8d %8.2f %8.2f %8.2f %8.2f\n',SPS(i,:));
% 输入深度、初始压力水头或土壤水分、根水分吸收的深度分布特征值（如果节点位于根区之外，设为0）
% 
end
fprintf(file_ID,'    0\n');
fclose(file_ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% write hydrus  %%%%%%%%%%%%%%

file_ID=fopen(strcat(file_path,'\','HYDRUS1D.DAT'),'wt');
fprintf(file_ID,';\n');
fprintf(file_ID,'[Main]\n');
fprintf(file_ID,'HYDRUS_Version=4\n');
fprintf(file_ID,'WaterFlow=1\n');
fprintf(file_ID,'SoluteTransport=0\n');
fprintf(file_ID,'Unsatchem=0\n');
fprintf(file_ID,'HP1=0\n');
fprintf(file_ID,'HeatTransport=0\n');
fprintf(file_ID,'EquilibriumAdsorption=1\n');
fprintf(file_ID,'MobileImmobile=0\n');
fprintf(file_ID,'RootWaterUptake=1\n');
fprintf(file_ID,'RootGrowth=1\n');
fprintf(file_ID,'MaterialNumbers=1\n');
fprintf(file_ID,'SubregionNumbers=1\n');
fprintf(file_ID,'SpaceUnit=cm\n');
fprintf(file_ID,'TimeUnit=hours\n');
fprintf(file_ID,'PrintTimes=1\n');
fprintf(file_ID,'NumberOfSolutes=0\n');
fprintf(file_ID,'InitialCondition=1\n');
fprintf(file_ID,';\n');
fprintf(file_ID,'[Profile]\n');
fprintf(file_ID,'NumberOfNodes=101\n');
fprintf(file_ID,'ProfileDepth=1.E+02\n');
% fprintf(file_ID,'ObservationNodes=0\n');
% fprintf(file_ID,'GridVisible=1\n');
% fprintf(file_ID,'SnapToGrid=1\n');
% fprintf(file_ID,'ProfileWidth=80\n');
% fprintf(file_ID,'LeftMargin=40\n');
% fprintf(file_ID,'GridOrgX=0\n');
% fprintf(file_ID,'GridOrgY=0\n');
% fprintf(file_ID,'GridDX=5.E+00\n');
% fprintf(file_ID,'GridDY=5.E+00\n');
fclose(file_ID);


fclose('all');
exe_str=[exe_path,' ',file_path];
dos(exe_str); % run HYDRUS-1D
quit cancel;
%==============================
% move outfile
%==============================
mkdir([file_path,'/T']);
mkdir([file_path,'/A']);
mkdir([file_path,'/P']);
mkdir([file_path,'/others']);
movefile([file_path,'/Obs_Node.out'],[file_path,'/T/Obs_Node.out']);
movefile([file_path,'/Run_Inf.out'],[file_path,'/T/Run_Inf.out']);
movefile([file_path,'/T_Level.out'],[file_path,'/T/T_Level.out']);
movefile([file_path,'/A_Level.out'],[file_path,'/A/A_Level.out']);
movefile([file_path,'/Meteo.out'],[file_path,'/A/Meteo.out']);
movefile([file_path,'/Balance.out'],[file_path,'/P/Balance.out']);
movefile([file_path,'/Nod_Inf.out'],[file_path,'/P/Nod_Inf.out']);
movefile([file_path,'/Nod_Inf_V.out'],[file_path,'/P/Nod_Inf_V.out']);
movefile([file_path,'/I_Check.out'],[file_path,'/others/I_Check.out']);
movefile([file_path,'/Profile.out'],[file_path,'/others/Profile.out']);
%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sd=100;
file_path='F:\HYDRUS-1D\work_dir\run_1\T\T_Level.out';
sim_data=importdata(file_path,' ',9);
sim_data=sim_data.data;

sim_sm=sim_data(:,17)./sd;
sim_evap=sim_data(:,4);
sim_tran=sim_data(:,5);
time=sim_data(:,1);
script=find(ismember(time,BC(:,1)));
sim_sm=sim_sm(script);
sim_evap=sim_evap(script);
sim_tran=sim_tran(script);

sim_et=sim_tran+sim_evap;
sim_et=sim_et*10*24*28.35;
sim_et(find(sim_evap <0 | sim_tran <0))=0;
script3=find(sim_et >1);
ET_data=[ET(script3),ET_balance(script3),sim_et(script3)];
ETEC_R2=corrcoef(ET(script3),sim_et(script3));
ETEC_R2=ETEC_R2(1,2)*ETEC_R2(1,2);
ETEC_bias=mean(sim_et(script3)-ET(script3));
ETEC_RMSE=sqrt(mean((sim_et(script3)-ET(script3)).^2));

ETEB_R2=corrcoef(ET_balance(script3),sim_et(script3));
ETEB_R2=ETEB_R2(1,2)*ETEB_R2(1,2);
ETEB_bias=mean(sim_et(script3)-ET_balance(script3));
ETEB_RMSE=sqrt(mean((sim_et(script3)-ET_balance(script3)).^2));

ET_tj(1,1)=ETEC_R2;ET_tj(1,2)=ETEC_bias;ET_tj(1,3)=ETEC_RMSE;
ET_tj(2,1)=ETEB_R2;ET_tj(2,2)=ETEB_bias;ET_tj(2,3)=ETEB_RMSE;

sm=0.1*sm(:,1)+0.15*sm(:,2)+0.25*sm(:,3)+0.5*sm(:,4);
script4=find(~isnan(sm));
sm=sm(script4);
sim_sm=sim_sm(script4);
rsm_data=[sm./100,sim_sm];

rsm_R2=corrcoef(sm./100,sim_sm);
rsm_R2=rsm_R2(1,2)*rsm_R2(1,2);
rsm_bias=mean(sim_sm-sm./100);
rsm_RMSE=sqrt(mean((sim_sm-sm./100).^2));
rsm_tj=[rsm_R2,rsm_bias,rsm_RMSE];

scatter(sm./100,sim_sm);
hold on
plot([0.1,0.5],[0.1,0.5],Color='red');
hold off
set(gca,'XLim',[0.1,0.5]);
set(gca,'YLim',[0.1,0.5]);
saveas(gcf,[num2str(year),'_rsm.jpg']);
close(gcf);
