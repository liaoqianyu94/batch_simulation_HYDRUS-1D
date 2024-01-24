function [BC,MC,albedo,sm,ET,ET_balance]=read_csv(obs_csv_path,simulation_strart_time,simulation_end_time)

% clc
% clear all
% 
% V_inform_path='F:\simulation_ameriflux\V_data.xlsx';
% site='US-Ne3';
% year=2007;
% [p_height,p_LAI,p_RD,year_span]=read_V(V_inform_path,site,year);
% 
% obs_csv_path='F:\simulation_ameriflux\AMF_US-Ne3_BASE-BADM_13-5\AMF_US-Ne3_BASE_HR_13-5.csv';
% simulation_strart_time=year_span(1)*10000;
% simulation_end_time=(year_span(2)+1)*10000;

obs=importdata(obs_csv_path);
obs_name=obs.colheaders;
obs_data=obs.data;

time=obs_data(:,1);
script=find(time >= simulation_strart_time & time < simulation_end_time);

Ta=obs_data(script,find(ismember(obs_name,'TA_PI_F_1_1_1')));

SW_IN=obs_data(script,find(ismember(obs_name,'SW_IN_PI_F_1_1_1')));
SW_OUT=obs_data(script,find(ismember(obs_name,'SW_OUT_PI_F_1_1_1'))); 
SW_OUT(find(SW_OUT==-9999))=0;

LW_IN=obs_data(script,find(ismember(obs_name,'LW_IN_PI_F_1_1_1')));
LW_OUT=obs_data(script,find(ismember(obs_name,'LW_OUT_PI_F_1_1_1')));
Net_R=SW_IN-SW_OUT+LW_IN-LW_OUT;

% script2=find(Net_R > 1000 | Net_R <-1000)
% for i=1:length(script2)
%     Net_R(script2(i))=Net_R(script2(i)-1)/2+Net_R(script2(i)+1)/2;
% end

albedo=SW_OUT./SW_IN;
albedo(find(albedo >1 | albedo <0))=0.23;
WS=obs_data(script,find(ismember(obs_name,'WS_1_2_1')));
script2=find(WS<0);
for i=1:length(script2)
    WS(script2(i))=WS(script2(i)-3)/2+WS(script2(i)+3)/2;
end

VPD=obs_data(script,find(ismember(obs_name,'RH_PI_F_1_1_1')));
Precipitation=obs_data(script,find(ismember(obs_name,'P_PI_F_2_2_1')));

%%%%%%%%%%%%%%%%% read sm  %%%%%%%%%%%%%%%%%%%%%%
sm_std=100;
sm_1=obs_data(script,find(ismember(obs_name,'SWC_PI_F_1_1_1')));
sm_2=obs_data(script,find(ismember(obs_name,'SWC_PI_F_2_1_1')));
sm_3=obs_data(script,find(ismember(obs_name,'SWC_PI_F_3_1_1')));

for i=1:length(sm_1)
L=[sm_1(i);sm_2(i);sm_3(i)]; 
L(find(L <0))=[];
if length(L)==3 & std(L)<sm_std
sm10(i,1)=mean(L);
else
sm10(i,1)=NaN;   
end
end

sm_1=obs_data(script,find(ismember(obs_name,'SWC_PI_F_1_2_1')));
sm_2=obs_data(script,find(ismember(obs_name,'SWC_PI_F_2_2_1')));
sm_3=obs_data(script,find(ismember(obs_name,'SWC_PI_F_3_2_1')));

for i=1:length(sm_1)
L=[sm_1(i);sm_2(i);sm_3(i)]; 
L(find(L <0))=[];
if length(L)==3 & std(L)<sm_std
sm25(i,1)=mean(L);
else
sm25(i,1)=NaN;   
end
end

sm_1=obs_data(script,find(ismember(obs_name,'SWC_PI_F_1_3_1')));
sm_2=obs_data(script,find(ismember(obs_name,'SWC_PI_F_2_3_1')));
sm_3=obs_data(script,find(ismember(obs_name,'SWC_PI_F_3_3_1')));

for i=1:length(sm_1)
L=[sm_1(i);sm_2(i);sm_3(i)]; 
L(find(L <0))=[];
if length(L)==3 & std(L)<sm_std
sm50(i,1)=mean(L);
else
sm50(i,1)=NaN;   
end
end

sm_1=obs_data(script,find(ismember(obs_name,'SWC_PI_F_1_4_1')));
sm_2=obs_data(script,find(ismember(obs_name,'SWC_PI_F_2_4_1')));
sm_3=obs_data(script,find(ismember(obs_name,'SWC_PI_F_3_4_1')));

for i=1:length(sm_1)
L=[sm_1(i);sm_2(i);sm_3(i)]; 
L(find(L <0))=[];
if length(L)==3 & std(L)<sm_std
sm100(i,1)=mean(L);
else
sm100(i,1)=NaN;   
end
end

sm=[sm10,sm25,sm50,sm100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
ET=obs_data(script,find(ismember(obs_name,'LE_PI_F_1_1_1')));
H=obs_data(script,find(ismember(obs_name,'H_PI_F_1_1_1')));
G=obs_data(script,find(ismember(obs_name,'G_PI_F_1_1_1')));
ET_balance=Net_R-H-G;

simulation_date=fix(time(script,1)/10000);
for i=1:length(simulation_date)
    Ta_max(i)=max(Ta(find(simulation_date==simulation_date(i))));
end

for i=1:length(simulation_date)
    Ta_min(i)=min(Ta(find(simulation_date==simulation_date(i))));
end

Ta_max=Ta_max';
Ta_min=Ta_min';

time_doy=ceil(juliandate(datetime(num2str(simulation_end_time),'InputFormat','yyyyMMddHHmm'))...
             -juliandate(datetime(num2str(simulation_strart_time),'InputFormat','yyyyMMddHHmm')));
time_num=time_doy*24;
if length(script) == time_num
    time_L=1:1:time_num;
    time_L=time_L';
elseif length(script) == time_num*2
    time_L=0.5:0.5:time_num;
    time_L=time_L';
end

rSoil=zeros(length(script),1);
rRoot=zeros(length(script),1);
hCritA=repelem(100000,length(script),1);
rB=zeros(length(script),1);
hB=zeros(length(script),1);
ht=zeros(length(script),1);

Precipitation=Precipitation*0.1;
BC=[time_L,Precipitation,rSoil,rRoot,hCritA,rB,hB,ht];


uniform=repelem(4,length(script),1);
WS=WS*86.4;
% SW_IN=SW_IN*0.0864;
Net_R=Net_R*0.0864;
VPD=VPD;
MC=[time_L,Net_R,Ta_max,Ta_min,VPD,WS,uniform];