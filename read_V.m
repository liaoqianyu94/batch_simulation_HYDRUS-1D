function [p_height,p_LAI,p_RD,year_span]=read_V(V_inform_path,site,year)

data_RD_max=120;
RD_start=5;
grass_height=50;
grass_RD=80;

if strcmp(site,'US-Ne1')
V_inform=xlsread(V_inform_path,'sheet1');
elseif strcmp(site,'US-Ne2')
V_inform=xlsread(V_inform_path,'sheet2');
elseif strcmp(site,'US-Ne3')
V_inform=xlsread(V_inform_path,'sheet3');
elseif strcmp(site,'US-Wkg')
V_inform=xlsread(V_inform_path,'sheet4');
end

    if strcmp(site,'US-Ne1') | strcmp(site,'US-Ne2') | strcmp(site,'US-Ne3')
    i=year-2002+1;
    N_start_year=V_inform(i,37);
    N_end_year=max(V_inform(:,2*i-1));
    year_span=[N_start_year;N_end_year];
   
   height_date=V_inform(1:250,2*i-1);
   script=find(height_date >0);
   height_date=height_date(script);
   height=V_inform(1:250,2*i);
   height=height(script);
   height=height-min(height)+0.1;
   height_span=length(height);
   
   emer_time=ceil(juliandate(datetime(num2str(min(height_date)),'InputFormat','yyyyMMdd'))...
             -juliandate(datetime(num2str(N_start_year),'InputFormat','yyyyMMdd')));
   emer_height=linspace(0,0,emer_time*24);
   
    V_growth_doy=datetime(num2str(height_date),'InputFormat','yyyyMMdd');
    V_growth_doy=day(V_growth_doy,'dayofyear');
    t= V_growth_doy-min(V_growth_doy);
    growth_height=interp1(t,height,1/24:1/24:max(t)+1);
    p_height=[emer_height';growth_height'];
    script=find(isnan(p_height));
    p_height(isnan(p_height))=p_height(script(1)-1);

    p_RD=-data_RD_max.*p_height./max(p_height);
    p_RD(find(p_RD<0))=p_RD(find(p_RD<0))-RD_start;

   LAI_date=V_inform(1:250,2*i+11);
   script=find(LAI_date >0);
   LAI_date=LAI_date(script);
   LAI=V_inform(1:250,2*i+12);
   LAI=LAI(script);
   
   emer_time=ceil(juliandate(datetime(num2str(min(LAI_date)),'InputFormat','yyyyMMdd'))...
             -juliandate(datetime(num2str(N_start_year),'InputFormat','yyyyMMdd')));
   emer_LAI=linspace(0,0,emer_time*24);

    V_growth_doy=datetime(num2str(LAI_date),'InputFormat','yyyyMMdd');
    V_growth_doy=day(V_growth_doy,'dayofyear');
    t= V_growth_doy-min(V_growth_doy);
    growth_LAI=interp1(t,LAI,1/24:1/24:max(t)+1);
    p_LAI=[emer_LAI';growth_LAI'];
    script=find(isnan(p_LAI));
    p_LAI(isnan(p_LAI))=p_LAI(script(1)-1);
  
    end

    if strcmp(site,'US-Wkg')
     i=year-2005+1
     LAI_date=V_inform(1:380,2*i+11);
     LAI_date=LAI_date(~isnan(LAI_date));
     LAI=V_inform(1:380,2*i+12);
     LAI=LAI(~isnan(LAI_date));

    V_growth_doy=datetime(num2str(LAI_date),'InputFormat','yyyyMMdd');
    V_growth_doy=day(V_growth_doy,'dayofyear');
    t= V_growth_doy-min(V_growth_doy);
    growth_LAI=interp1(t,LAI,1/24:1/24:max(t)+1);
    p_LAI=growth_LAI';
    
    p_height(1:length(p_LAI),1)=grass_height;
    p_RD(1:length(p_LAI),1)=-grass_RD;
       
    end
