function DisplayCloudTopPhaseMovie(titlestr,iframe)
% Display a single movie frame of the cloud to phase data
 
% Written By: Stephen Forczyk
% Created: Mar 30,2022
% Revised: --------
% Classification: Unclassified

global BandDataS MetaDataS GoesWaveBand;
global RadS DQFS tS yS xS tBS goesImagerS yImgS yImgBS;
global xImgS xImgBS SatDataS GeoSpatialS;
global PhaseS AlgoS ErrorS VersionContainerS;
global Algo1S ProcessParamS;
global SZAS SZABS Error1S CloudPixelsS;
global GOESFileName itype;
global UnSatPixS StarIDS;
global RptGenPresent iCreatePDFReport pdffilename rpt chapter ifirstCall;
global JpegCounter JpegFileList;
global DQFHdr DQFTable MiscHdr MiscTable RqmtsHdr RqmtsTable;
global PhaseHdr PhaseTable;
global numncfiles;
global MonthDayStr MonthDayYearStr;
global DayMonthNonLeapYear DayMonthLeapYear CalendarFileName;
global westEdge eastEdge northEdge southEdge;
global MapFormFactor;
global RasterLats RasterLons;
global FlashDataByFrame GroupPosByFrame;
global GroupMaxFlashDurations FileList GroupLongestFlashData ifirstCallGLMLightning;
global GMTCTPFracS GMTCTPFracE;
global iaddLightning FOVLOC MOVLOC;
global TempMovieName vTemp numMovieFrames nowMovieFrame;
global GridMatFileName agrid;
global PhaseHistory;

global widd2 lend2;
global initialtimestr igrid ijpeg ilog imovie;
global vert1 hor1 widd lend;
global vert2 hor2 machine;
global chart_time;
global Fz1 Fz2 fid;
global idirector mov izoom iwindow;
global matpath GOES16path;
global jpegpath ;
global smhrpath excelpath ascpath;
global ipowerpoint PowerPointFile scaling stretching padding;
global ichartnum;
% additional paths needed for mapping
global matpath1 mappath;
global canadapath stateshapepath topopath;
global trajpath militarypath gridpath;
global figpath screencapturepath;
global shapepath2 countrypath countryshapepath usstateboundariespath;
global GOES16Band1path GOES16Band2path GOES16Band3path GOES16Band4path;
global GOES16Band5path GOES16Band6path GOES16Band7path GOES16Band8path
global GOES16Band9path GOES16Band10path GOES16Band11path GOES16Band12path;
global GOES16Band13path GOES16Band14path GOES16Band15path GOES16Band16path;
global GOES16CloudTopHeightpath;

if((iCreatePDFReport==1) && (RptGenPresent==1))
    import mlreportgen.dom.*;
    import mlreportgen.report.*;
end
if(iframe==1)
    fprintf(fid,'\n');
    if(itype==1)
        fprintf(fid,'%s\n','------- Plot Cloud Top Phase Data For Full Disk ------');
    elseif(itype==2)
        fprintf(fid,'%s\n','------- Plot Cloud Top Phase Data For Conus ------');
    elseif(itype==3)
        fprintf(fid,'%s\n','------- Plot Cloud Top Phase Data Meso1 Scale ------'); 
    elseif(itype==4)
        fprintf(fid,'%s\n','------- Plot Cloud Top Phase Data Meso2 Scale ------'); 
    else

    end
end
eval(['cd ' gridpath(1:length(gridpath)-1)]);
ishiftnow=FOVLOC{1+iframe,13};
frameint=10;
if(numncfiles>500)
    frameint=40;
end
rm=mod(iframe-1,frameint);
if(iframe==1)
    agrid=-1;
end
if((iframe==1) && (itype==1))
    GridMatFileName='FDCloudTopPhaseBoundaries.mat';
    agrid=exist(GridMatFileName,'file'); 
    waitstr='Calculating Full Disk Geographic Raster';
elseif((iframe==1) && (itype==2))
    GridMatFileName='ConusCloudTopPhaseBoundaries.mat'; 
    agrid=exist(GridMatFileName,'file'); 
    waitstr='Calculating Conus Geographic Raster';
elseif((iframe==1) && (itype==3))
    GridMatFileName='Meso1CloudTopPhaseBoundaries.mat';
    agrid=exist(GridMatFileName,'file');
    waitstr='Calculating Meso1 Geographic Raster';
elseif((iframe==1) && (itype==4))
    GridMatFileName='Meso2CloudTopPhaseBoundaries.mat'; 
    agrid=exist(GridMatFileName,'file'); 
    waitstr='Calculating Meso2 Geographic Raster';
end
if((ishiftnow==1) && (itype==1))
    GridMatFileName='FDCloudTopPhaseBoundaries.mat';
    waitstr='Re-Calculating Full Disk Geographic Raster';  
    agrid=0;
elseif((ishiftnow==1) && (itype==2))
    GridMatFileName='ConusCloudTopPhaseBoundaries.mat';  
    waitstr='Re-Calculating Conus Geographic Raster'; 
    agrid=0;
elseif((ishiftnow==1) && (itype==3))
    GridMatFileName='Meso1CloudTopPhaseBoundaries.mat';
    waitstr='Re-Calculating Meso1 Geographic Raster'; 
    agrid=0;
elseif((ishiftnow==1) && (itype==4))
    GridMatFileName='Meso2CloudTopPhaseBoundaries.mat';
    waitstr='Re-Calculating Meso2 Geographic Raster';   
    agrid=0;
end
if(agrid~=2)
    numrasterlon=length(xS.values);
    numrasterlat=length(yS.values);
    xx=xS.values;
    yy=yS.values;
    RasterLats=zeros(numrasterlat,numrasterlon);
    RasterLons=RasterLats;
    h=waitbar(0,waitstr);
    tic;
    for jj=1:numrasterlat
        for kk=1:numrasterlon
            [GeodLat,GeodLon] = CalculateGeodCoordFromXYPosRev1(xx(kk,1),yy(jj,1));
            RasterLats(jj,kk)=real(GeodLat);% This okay even in reverse order they will be transposed later
            RasterLons(jj,kk)=real(GeodLon);% This okay even in reverse order they will be transposed later
        end
        waitbar(jj/numrasterlat);
    end
    close(h);
    elapsed_time=toc;
    dispstr=strcat('Calculating the Raster took-',num2str(elapsed_time,5),'-sec');
    disp(dispstr)
    actionstr='save';
    varstr='RasterLats RasterLons';
    qualstr='-v7.3';
    [cmdString]=MyStrcatV73(actionstr,GridMatFileName,varstr,qualstr);
    eval(cmdString)
    dispstr=strcat('Saved Raster Data To Matlab File-',GridMatFileName);
    disp(dispstr);
    agrid=2;
end
if(agrid~=2)
    gridstr1=strcat('Loaded Map Grid From File=',GridMatFileName);
    fprintf(fid,'%s\n',gridstr1);
    load(GridMatFileName,'RasterLats','RasterLons');
    dispstr=strcat('Loaded Grid Raster File-',GridMatFileName);
    disp(dispstr)
end
Phase=PhaseS.values;
measunits=PhaseS.units;
scalebarlen=1000;
scalebarloc='nw'; % location of scale (sw=southwest)
[nrows,ncols]=size(Phase);
minRasterLat=min(min(RasterLats));
minRasterLon=min(min(RasterLons));
maxRasterLat=max(max(RasterLats));
maxRasterLon=max(max(RasterLons));
if(iframe==1)
    gridstr2='****Grid Limits Follow*****';
    fprintf(fid,'%s\n',gridstr2);
    gridstr1=strcat('Will Load Map Grid From File=',GridMatFileName);
    fprintf(fid,'%s\n',gridstr1);

    gridstr3=strcat('Min Raster Lat=',num2str(minRasterLat,6),'-Max RasterLat=',num2str(maxRasterLat,6));
    fprintf(fid,'%s\n',gridstr3);
    gridstr4=strcat('Min Raster Lon=',num2str(minRasterLon,6),'-Max RasterLon=',num2str(maxRasterLon,6));
    fprintf(fid,'%s\n',gridstr4);
    gridstr6='****Grid Limits End*****';
    fprintf(fid,'%s\n',gridstr6);
end
% Get some stats on the Phase
Phase1D=reshape(Phase,nrows*ncols,1);
Phase1DS=sort(Phase1D);
maxPhase=max(Phase1D);
minPhase=min(Phase1D);
nvals=nrows*ncols;
NaNPhase=isnan(Phase1D);
numnan=sum(NaNPhase);
nanfrac=numnan/(nrows*ncols);
validfrac=1-nanfrac;
% Fetch the map limits
westEdge=double(GeoSpatialS.geospatial_westbound_longitude);
eastEdge=double(GeoSpatialS.geospatial_eastbound_longitude);
southEdge=double(GeoSpatialS.geospatial_southbound_latitude);
northEdge=double(GeoSpatialS.geospatial_northbound_latitude);
if(iframe==1)
    maplimitstr1='****Map Limits Follow*****';
    fprintf(fid,'%s\n',maplimitstr1);
    maplimitstr2=strcat('WestEdge=',num2str(westEdge,7),'-EastEdge=',num2str(eastEdge));
    fprintf(fid,'%s\n',maplimitstr2);
    maplimitstr3=strcat('SouthEdge=',num2str(southEdge,7),'-NorthEdge=',num2str(northEdge));
    fprintf(fid,'%s\n',maplimitstr3);
    maplimitstr4='****Map Limits End*****';
    fprintf(fid,'%s\n',maplimitstr4);
end

% Get the scan start date and time data based on the file name
[YearS,DayS,HourS,MinuteS,SecondS] = GetGOESDateTimeString(GOESFileName,1);
% Get the scan end date and time data based on the file name
[YearE,DayE,HourE,MinuteE,SecondE] = GetGOESDateTimeString(GOESFileName,2);
result = is_leap_year(YearS);
if(result==0)
    MonthDayStr=char(DayMonthNonLeapYear{DayS,1});
else
    MonthDayStr=char(DayMonthLeapYear{DayS,1});
end
MonthDayYearStr=strcat(MonthDayStr,'-',num2str(YearS));
result = is_leap_year(YearS);
% Get the scan start date and time data based on the file name
GMTCTPIS=60*MinuteS+SecondS;
GMTCTPFracS=HourS+GMTCTPIS/3600;
GMTCTPIE=60*MinuteE+SecondE;
GMTCTPFracE=HourE+GMTCTPIE/3600;
if(iframe==1)
    datestr10='***** Scan Time Definition Follows *****';
    datestr14='***** Scan Time Definition End *****';
    fprintf(fid,'%s\n',datestr10);
    datestr11=strcat('Year=',num2str(YearS,4),'-Day=',num2str(DayS,3),'-Hour=',num2str(HourS,2),...
        '-Minute=',num2str(MinuteS,2),'-Second=',num2str(SecondS,5));
    fprintf(fid,'%s\n',datestr11);
    datestr12=strcat('Calender Month and Day=',MonthDayStr);
    fprintf(fid,'%s\n',datestr12);
    fprintf(fid,'%s\n',datestr14);
end
% Get a count of the phase by pixels
Phase0=floor(nvals*PhaseS.percent_clear_sky);
Phase0pct=PhaseS.percent_clear_sky;
Phase0fracstr=num2str(Phase0pct,5);
Phase1=floor(nvals*PhaseS.percent_liquid_water);
Phase1pct=PhaseS.percent_liquid_water;
Phase1fracstr=num2str(Phase1pct,5);
Phase2=floor(nvals*PhaseS.percent_super_cooled_liquid_water);
Phase2pct=PhaseS.percent_super_cooled_liquid_water;
Phase2fracstr=num2str(Phase2pct,5);
Phase3=floor(nvals*PhaseS.percent_mixed_phase);
Phase3pct=PhaseS.percent_mixed_phase;
Phase3fracstr=num2str(Phase3pct,5);
Phase4=floor(nvals*PhaseS.percent_ice);
Phase4pct=PhaseS.percent_ice;
Phase4fracstr=num2str(Phase4pct,5);
Phase5=floor(nvals*PhaseS.percent_unknown);
Phase5pct=PhaseS.percent_unknown;
Phase5fracstr=num2str(Phase5pct,5);
phasestr0=strcat('Phase 0-% clear sky=',num2str(Phase0pct,5),'-pixels=',num2str(Phase0));
phasestr1=strcat('Phase 1-% liquid water=',num2str(Phase1pct,5),'-pixels=',num2str(Phase1));
phasestr2=strcat('Phase 2-% super cooled liquid water=',num2str(Phase2pct,5),'-pixels=',num2str(Phase2));
phasestr3=strcat('Phase 3-% mixed phase=',num2str(Phase3pct,5),'-pixels=',num2str(Phase3));
phasestr4=strcat('Phase 4-% ice=',num2str(Phase4pct,5),'-pixels=',num2str(Phase4));
phasestr5=strcat('Phase 5-% unknown phase=',num2str(Phase5pct,5),'-pixels=',num2str(Phase5));
phasestr6=strcat('nanfrac=',num2str(nanfrac,5),'-validfrac=',num2str(validfrac,5));
if(iframe==1)
    fprintf(fid,'%s\n',phasestr0);
    fprintf(fid,'%s\n',phasestr1);
    fprintf(fid,'%s\n',phasestr2);
    fprintf(fid,'%s\n',phasestr3);
    fprintf(fid,'%s\n',phasestr4);
    fprintf(fid,'%s\n',phasestr5);
    fprintf(fid,'%s\n',phasestr6);
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','    Frame    Time   Clear Sky Frac   Liquid Water  SCool Water    Mixed Phase   Ice Phase Unknown Phase');
end
fprintf(fid,'     %d  %10.4f   %10.5f  ',iframe,GMTCTPFracS,Phase0pct);
fprintf(fid,'   %10.5f   %10.5f   ',Phase1pct,Phase2pct);
fprintf(fid,'   %10.5f   %10.5f   ',Phase3pct,Phase4pct);
fprintf(fid,'%10.5f\n',Phase5pct);
%% Add data to the PhaseHistory Cell
PhaseHistory{1+iframe,1}=iframe;
PhaseHistory{1+iframe,2}=GOESFileName;
PhaseHistory{1+iframe,3}=GMTCTPFracS;
PhaseHistory{1+iframe,4}=GMTCTPFracE;
PhaseHistory{1+iframe,5}=nanfrac;
PhaseHistory{1+iframe,6}=Phase0pct;
PhaseHistory{1+iframe,7}=Phase1pct;
PhaseHistory{1+iframe,8}=Phase2pct;
PhaseHistory{1+iframe,9}=Phase3pct;
PhaseHistory{1+iframe,10}=Phase4pct;
PhaseHistory{1+iframe,11}=Phase5pct;
%% Now if the add lightning option is selected see if there is a file that matches this time period that has the data
% on the longest lasting Group by Time
igotbigdata=0;
if(iaddLightning==1)
    [nrowsldp1,ncols]=size(FlashDataByFrame);
    nrowsld=nrowsldp1-1;
    igo2=1;
    ik=0;
    while(igo2>0)
        ik=ik+1;
        TS=FlashDataByFrame{1+ik,9};
        TE=FlashDataByFrame{1+ik,10};
        if(itype<3)
            if((TS>=GMTCTPFracS) && (TE<=GMTCTPFracE))
                igo2=0;
                GroupBigLats=GroupLongestFlashData{1+ik,5};
                GroupBigLons=GroupLongestFlashData{1+ik,6};
                igotbigdata=1;
                iksave1=ik;
            end
            if(ik>=nrowsld)
                igo2=0;
            end
        elseif(itype>2)
           if((TS>=GMTCTPFracS) && (TE>=GMTCTPFracE))
                igo2=0;
                GroupBigLats=GroupLongestFlashData{1+ik,5};
                GroupBigLons=GroupLongestFlashData{1+ik,6};
                igotbigdata=1;
                iksave1=ik;
            end
            if(ik>=nrowsld)
                igo2=0;
            end
        end
    end
end
%% Now retrieve the data for each group by frame-not just the longest lasting groups
igotgroupdata=0;
if(iaddLightning==1)
    [nrowsldp1,ncols]=size(GroupPosByFrame);
    nrowsld=nrowsldp1-1;
    igo2=1;
    ik=0;
    iksave=0;
    while(igo2>0)
        ik=ik+1;
        TS=GroupPosByFrame{1+ik,7};
        TE=GroupPosByFrame{1+ik,8};
        if((TS>=GMTCTPFracS) && (TE>=GMTCTPFracE))
            igo2=0;
            GroupLats=GroupPosByFrame{1+ik,5};
            GroupLons=GroupPosByFrame{1+ik,6};
            igotgroupdata=1;
            iksave2=ik;
        end
        if(ik>=nrowsld)
            igo2=0;
        end
    end
end
% Set up the map axis
if(itype==1)
    axesm ('globe','Frame','on','Grid','on','meridianlabel','off','parallellabel','on',...
        'plabellocation',[-60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60],'mlabellocation',[]);
elseif(itype==2)
    axesm ('pcarree','Frame','on','Grid','on','MapLatLimit',[southEdge northEdge],...
     'MapLonLimit',[westEdge eastEdge],'meridianlabel','on','parallellabel','on','plabellocation',10,'mlabellocation',20,...
     'MLabelParallel','south','Frame','on','FontColor','b'); 
     gridm('GLineStyle','-','Gcolor',[.1 .1 .1],'Galtitude',.002,'MLineLocation',10,...
    'PLineLocation',10)
elseif((itype==3)||(itype==4))
    axesm ('pcarree','Frame','on','Grid','on','MapLatLimit',[southEdge northEdge],...
     'MapLonLimit',[westEdge eastEdge],'meridianlabel','on','parallellabel','on','plabellocation',2,'mlabellocation',2,...
     'MLabelParallel','south','Frame','on','FontColor','b');
    gridm('GLineStyle','-','Gcolor',[.1 .1 .1],'Galtitude',.002,'MLineLocation',2,...
    'PLineLocation',2)
end
set(gcf,'MenuBar','none');
set(gcf,'Position',[hor1 vert1 widd lend])
geoshow(RasterLats,RasterLons,Phase,'DisplayType','texturemap');
zlimits=[0 5];
demcmap(zlimits);
hold on
% load the country borders and plot them
eval(['cd ' mappath(1:length(mappath)-1)]);
load('USAHiResBoundaries.mat','USALat','USALon');
plotm(USALat,USALon,'r');
load('CanadaBoundaries.mat','CanadaLat','CanadaLon');
plotm(CanadaLat,CanadaLon,'r');
load('MexicoBoundaries.mat','MexicoLat','MexicoLon');
plotm(MexicoLat,MexicoLon,'r');
load('CubaBoundaries.mat','CubaLat','CubaLon');
plotm(CubaLat,CubaLon,'r');
load('DominicanRepublicBoundaries.mat','DRLat','DRLon');
plotm(DRLat,DRLon,'r');
load('HaitiBoundaries.mat','HaitiLat','HaitiLon');
plotm(HaitiLat,HaitiLon,'r');
load('BelizeBoundaries.mat','BelizeLat','BelizeLon');
plotm(BelizeLat,BelizeLon,'r');
load('GautemalaBoundaries.mat','GautemalaLat','GautemalaLon');
plotm(GautemalaLat,GautemalaLon,'r')
load('JamaicaBoundaries.mat','JamaicaLat','JamaicaLon');
plotm(JamaicaLat,JamaicaLon,'r');
load('NicaraguaBoundaries.mat','NicaraguaLat','NicaraguaLon');
plotm(NicaraguaLat,NicaraguaLon,'r')
load('HondurasBoundaries.mat','HondurasLat','HondurasLon');
plotm(HondurasLat,HondurasLon,'r')
load('ElSalvadorBoundaries.mat','ElSalvadorLat','ElSalvadorLon');
plotm(ElSalvadorLat,ElSalvadorLon,'r');
load('PanamaBoundaries.mat','PanamaLat','PanamaLon');
plotm(PanamaLat,PanamaLon,'r');
load('ColumbiaBoundaries.mat','ColumbiaLat','ColumbiaLon');
plotm(ColumbiaLat,ColumbiaLon,'r');
load('VenezuelaBoundaries.mat','VenezuelaLat','VenezuelaLon');
plotm(VenezuelaLat,VenezuelaLon,'r')
load('PeruBoundaries.mat','PeruLat','PeruLon');
plotm(PeruLat,PeruLon,'r');
load('EcuadorBoundaries.mat','EcuadorLat','EcuadorLon');
plotm(EcuadorLat,EcuadorLon,'r')
load('BrazilBoundaries.mat','BrazilLat','BrazilLon');
plotm(BrazilLat,BrazilLon,'r');
load('BoliviaBoundaries.mat','BoliviaLat','BoliviaLon');
plotm(BoliviaLat,BoliviaLon,'r')
load('ChileBoundaries.mat','ChileLat','ChileLon');
plotm(ChileLat,ChileLon,'r');
load('ArgentinaBoundaries.mat','ArgentinaLat','ArgentinaLon');
plotm(ArgentinaLat,ArgentinaLon,'r');
load('UruguayBoundaries.mat','UruguayLat','UruguayLon');
plotm(UruguayLat,UruguayLon,'r');
load('CostaRicaBoundaries.mat','CostaRicaLat','CostaRicaLon');
plotm(CostaRicaLat,CostaRicaLon,'r');
load('FrenchGuianaBoundaries.mat','FrenchGuianaLat','FrenchGuianaLon');
plotm(FrenchGuianaLat,FrenchGuianaLon,'r');
load('GuyanaBoundaries.mat','GuyanaLat','GuyanaLon');
plotm(GuyanaLat,GuyanaLon,'r');
load('SurinameBoundaries.mat','SurinameLat','SurinameLon');
plotm(SurinameLat,SurinameLon,'r');
load('PuertoRico_Boundaries.mat','PRLat','PRLon');
plotm(PRLat,PRLon,'r');
load('Turks_Boundaries.mat','TurksLat','TurksLon');
plotm(TurksLat,TurksLon,'r')
load('VirginIslandsBoundaries.mat','VirginIslandsLat','VirginIslandsLon');
plotm(VirginIslandsLat,VirginIslandsLon,'r')
load('BahamasBoundaries.mat','BahamasLat','BahamasLon');
plotm(BahamasLat,BahamasLon,'r')
load('BarbadosBoundaries.mat','BarbadosLat','BarbadosLon');
plotm(BarbadosLat,BarbadosLon,'r')
load('BritishVirginIslandsBoundaries.mat','BVALat','BVALon');
plotm(BVALat,BVALon,'r')
load('CaymanBoundaries.mat','CaymanLat','CaymanLon');
plotm(CaymanLat,CaymanLon,'r')
load('DominicaBoundaries.mat','DominicaLat','DominicaLon');
plotm(DominicaLat,DominicaLon,'r')
load('GrenadaBoundaries.mat','GrenadaLat','GrenadaLon');
plotm(GrenadaLat,GrenadaLon,'r')
load('GuadeloupeBoundaries.mat','GuadeloupeLat','GuadeloupeLon');
plotm(GuadeloupeLat,GuadeloupeLon,'r')
load('MartiniqueBoundaries.mat','MartiniqueLat','MartiniqueLon');
plotm(MartiniqueLat,MartiniqueLon,'r')
load('SaintHelenaBoundaries.mat','SaintHelenaLat','SaintHelenaLon');
plotm(SaintHelenaLat,SaintHelenaLon,'r')
load('SaintKittsBoundaries.mat','SaintKittsLat','SaintKittsLon');
plotm(SaintKittsLat,SaintKittsLon,'r')
load('SaintPierreBoundaries.mat','SaintPierreLat','SaintPierreLon');
plotm(SaintPierreLat,SaintPierreLon,'r')
load('SaintVincentBoundaries.mat','SaintVincentLat','SaintVincentLon');
plotm(SaintVincentLat,SaintVincentLon,'r')
if(itype==1)
    view(0,0);
end
title(titlestr)
tightmap;
if((itype==3)|| (itype==4))
    scalebar('length',scalebarlen,'units','km','color','r','location',scalebarloc,'fontangle','bold',...
        'RulerStyle','patches','orientation','vertical');
end
hc=colorbar;
%% Now plot in each Group Flash Location
if((iaddLightning==1)  && (igotgroupdata>0))
        plotm(GroupLats,GroupLons,'ro',...
        'MarkerSize',4,'MarkerEdgeColor',[1 1 1],...
        'MarkerFaceColor',[1 0.5 0]);
end
%% Now plot the longest duration event members
if((iaddLightning==1) && (igotbigdata>0))
        plotm(GroupBigLats,GroupBigLons,'rd',...
        'MarkerSize',6,'MarkerEdgeColor',[1 1 1],...
        'MarkerFaceColor',[0 1 1]);
end
hold off
newaxesh=axes('Position',[0 0 1 1]);
set(newaxesh,'XLim',[0 1],'YLim',[0 1]);
if(itype==1)
    Pos=get(hc,'Position');
    PosyL=.80*Pos(1,2);
    PosyU=1.03*(Pos(1,2)+Pos(1,4));
    tx1=.13;
    ty1=PosyL;
    tx2=tx1;
    ty2=PosyL;
    txtstr2=strcat('Phase 0- clear sky=',num2str(Phase0pct,5),'-Phase 1-Liquid Water=',num2str(Phase1pct,5));
    txt2=text(tx2,ty2,txtstr2,'FontWeight','bold','FontSize',12);
    tx2a=tx2;
    ty2a=ty2-.03;
    txtstr2a=strcat('Phase 2- SCool water=',num2str(Phase2pct,5),'-Phase 3-Mixed=',num2str(Phase3pct,5));
    txt2a=text(tx2a,ty2a,txtstr2a,'FontWeight','bold','FontSize',12);
    tx2b=tx2a;
    ty2b=ty2a-.03;
    txtstr2b=strcat('Phase 4- Ice=',num2str(Phase4pct,5),'-Phase 5-Unknown=',num2str(Phase5pct,5));
    txt2b=text(tx2b,ty2b,txtstr2b,'FontWeight','bold','FontSize',12);
    tx3=Pos(1,1);
    ty3=1.01*PosyU;
    txtstr3='Phase';
    txt3=text(tx3,ty3,txtstr3,'FontWeight','bold','FontSize',10);
    tx4=tx1;
    ty4=ty2-.04;
elseif(itype==2)
    Pos=get(hc,'Position');
    PosyL=.80*Pos(1,2);
    PosyU=1.03*(Pos(1,2)+Pos(1,4));
    tx1=.13;
    ty1=PosyL;
    tx2=tx1;
    ty2=PosyL;
    txtstr2=strcat('Phase 0- clear sky=',num2str(Phase0pct,5),'-Phase 1-Liquid Water=',num2str(Phase1pct,5));
    txt2=text(tx2,ty2,txtstr2,'FontWeight','bold','FontSize',12);
    tx2a=tx2;
    ty2a=ty2-.03;
    txtstr2a=strcat('Phase 2- SCool water=',num2str(Phase2pct,5),'-Phase 3-Mixed=',num2str(Phase3pct,5));
    txt2a=text(tx2a,ty2a,txtstr2a,'FontWeight','bold','FontSize',12);
    tx2b=tx2a;
    ty2b=ty2a-.03;
    txtstr2b=strcat('Phase 4- Ice=',num2str(Phase4pct,5),'-Phase 5-Unknown=',num2str(Phase5pct,5));
    txt2b=text(tx2b,ty2b,txtstr2b,'FontWeight','bold','FontSize',12);
    tx3=Pos(1,1);
    ty3=1.01*PosyU;
    txtstr3='Phase';
    txt3=text(tx3,ty3,txtstr3,'FontWeight','bold','FontSize',10);
    tx4=tx2b;
    ty4=ty2b-.03;
    txt4=text(tx4,ty4,phasestr6,'FontWeight','bold','FontSize',12);
elseif((itype==3)|| (itype==4))
    Pos=get(hc,'Position');
    PosyL=.80*Pos(1,2);
    PosyU=1.03*(Pos(1,2)+Pos(1,4));
    tx1=.13;
    ty1=PosyL;
    tx2=tx1;
    ty2=PosyL-.03;
    txtstr2=strcat('Phase 0- clear sky=',num2str(Phase0pct,5),'-Phase 1-Liquid Water=',num2str(Phase1pct,5));
    txt2=text(tx2,ty2,txtstr2,'FontWeight','bold','FontSize',10);
    tx2a=tx2;
    ty2a=ty2-.02;
    txtstr2a=strcat('Phase 2- SCool water=',num2str(Phase2pct,5),'-Phase 3-Mixed=',num2str(Phase3pct,5));
    txt2a=text(tx2a,ty2a,txtstr2a,'FontWeight','bold','FontSize',10);
    tx2b=tx2a;
    ty2b=ty2a-.02;
    txtstr2b=strcat('Phase 4- Ice=',num2str(Phase4pct,5),'-Phase 5-Unknown=',num2str(Phase5pct,5));
    txt2b=text(tx2b,ty2b,txtstr2b,'FontWeight','bold','FontSize',10);
    tx3=Pos(1,1);
    ty3=1.01*PosyU;
    txtstr3='Phase';
    txt3=text(tx3,ty3,txtstr3,'FontWeight','bold','FontSize',10);

end
set(newaxesh,'Visible','Off');
% Save this chart
if(rm==1)
    figstr=strcat(titlestr,'.jpg');
    eval(['cd ' jpegpath(1:length(jpegpath)-1)]);
    actionstr='print';
    typestr='-djpeg';
    [cmdString]=MyStrcat2(actionstr,typestr,figstr);
    eval(cmdString);
%     savestr=strcat('Saved Cloud Top Phase data plot to file-',figstr);
%     fprintf(fid,'%s\n',savestr);
end
pause(chart_time);
% JpegCounter=JpegCounter+1;
% JpegFileList{JpegCounter+1,1}=figstr;
% JpegFileList{JpegCounter+1,2}=jpegpath;
frame=getframe(gcf);
writeVideo(vTemp,frame);
close('all');
if((iCreatePDFReport==1) && (RptGenPresent==1) && (rm==1))
    ifirstCall=ifirstCall+1;
    [ibreak]=strfind(GOESFileName,'_e');
    is=1;
    ie=ibreak(1)-1;
    GOESShortFileName=GOESFileName(is:ie);
    headingstr1=strcat('Cloud Top Phase-',GOESShortFileName);
    if(ifirstCall==1)
        chapter = Chapter("Title",headingstr1);
    end
    sectionstr=strcat('CloudPhase-',MapFormFactor);
    add(chapter,Section(sectionstr));
    imdata = imread(figstr);
    [nhigh,nwid,~]=size(imdata);
    image = mlreportgen.report.FormalImage();
    image.Image = which(figstr);% 
    pdftxtstr=strcat('Cloud Top Phase-For File-',GOESShortFileName);
    pdftext = Text(pdftxtstr);
    pdftext.Color = 'red';
    image.Caption = pdftext;
    nhighs=floor(nhigh/2.5);
    nwids=floor(nwid/2.5);
    heightstr=strcat(num2str(nhighs),'px');
    widthstr=strcat(num2str(nwids),'px');
    image.Height = heightstr;
    image.Width = widthstr;
    image.ScaleToFit=0;
    add(chapter,image);
% Now add some text -start by decribing the Cloud Top Phase Data
    parastr11=strcat('The Data for this chart was taken from file-',GOESShortFileName,'.');
    parastr12='Plotted on the chart is the Cloud Top Phase data as seen by the ABI sensor.';
    parastr13=strcat('This Map uses data taken in the-',MapFormFactor,'-format.');
    parastr14='Note that the phase data is reported by pixel with integer values ranging from 0 to 5.';
    parastr15='The detailed break out of the fraction of pixels in each of the 6 phases are shown in the next table.';
    parastr16='The algorithm that calculates the Cloud Top Phase uses data from Channels 10,11,14,and 15 and operates day or night.';
    if(iaddLightning==0)
        parastr17='Since the iaddLightning flag was set to zero lightning data was not added to the figure.';
    else
        parastr17='Since the iaddLightning flag was set to 1 lightning data was  added to the figure using colored circles and diamonds.';
    end
    parastr19=strcat(parastr11,parastr12,parastr13,parastr14,parastr15,parastr16,parastr17);
    p1 = Paragraph(parastr19);
    p1.Style = {OuterMargin("0pt", "0pt","10pt","10pt")};
    add(chapter,p1);
% Now build a phase Table to hold the results
    PhaseHdr={'Phase Value' 'Meaning'  'Pixel Frac'};
    PhaseTable=cell(6,3);
    PhaseTable{1,1}='0';
    PhaseTable{1,2}='Clear Sky';
    PhaseTable{1,3}=Phase0fracstr;
    PhaseTable{2,1}='1';
    PhaseTable{2,2}='Liquid Water';
    PhaseTable{2,3}=Phase1fracstr;
    PhaseTable{3,1}='2';
    PhaseTable{3,2}='Super Coooled Liquid Water';
    PhaseTable{3,3}=Phase2fracstr;
    PhaseTable{4,1}='3';
    PhaseTable{4,2}='Mixed Phase';
    PhaseTable{4,3}=Phase3fracstr;
    PhaseTable{5,1}='4';
    PhaseTable{5,2}='Ice';
    PhaseTable{5,3}=Phase4fracstr;
    PhaseTable{6,1}='5';
    PhaseTable{6,2}='Unknown';
    PhaseTable{6,3}=Phase5fracstr;
    T2=[PhaseHdr;PhaseTable];
    tb21=Table(T2);
    tb21.Style = [tb21.Style {Border('solid','black','3px')}];
    tb21.HAlign='center';
    tb21.TableEntriesHAlign = 'center';
    tb21.ColSep = 'single';
    tb21.RowSep = 'single';
    r = row(tb21,1);
    r.Style = [r.Style {Color('red'),Bold(true)}];
    bt2 = BaseTable(tb21);
    tabletitle = Text('Cloud Top Phase Table');
    tabletitle.Bold = false;
    bt2.Title = tabletitle;
    bt2.TableWidth="7in";
    add(chapter,bt2);
% Now build a small data quality table
    dqfflag1=100*DQFS.percent_overall_good_quality_qf;
    dqfflagstr1=num2str(dqfflag1,6);
    dqfflag2=100*DQFS.percent_overall_degraded_quality_qf;
    dqfflagstr2=num2str(dqfflag2,6);
    dqfflag3=100*DQFS.percent_good_quality_L1b_data_qf;
    dqfflagstr3=num2str(dqfflag3,6);
    dqfflag4=100*DQFS.percent_degraded_quality_L1b_data_qf;
    dqfflagstr4=num2str(dqfflag4,6);
    dqfflag5=100*DQFS.percent_good_quality_beta_ratio_qf;
    dqfflagstr5=num2str(dqfflag5,6);
    dqfflag6=100*DQFS.percent_degraded_quality_beta_ratio_qf;
    dqfflagstr6=num2str(dqfflag6,6);
    dqfflag7=100*DQFS.pct_ice_cloud_determ_based_on_strong_rad_signal_qf;
    dqfflagstr7=num2str(dqfflag7,6);
    dqfflag8=100*DQFS.pct_ice_cloud_determ_based_on_weak_rad_signal_qf;
    dqfflagstr8=num2str(dqfflag8,6);
    dqfflag9=100*DQFS.percent_good_quality_surface_emissivity_qf;
    dqfflagstr9=num2str(dqfflag9,6);
    dqfflag10=100*DQFS.percent_degraded_quality_surface_emissivity_qf;
    dqfflagstr10=num2str(dqfflag10,6);
    dqfflag11=100*DQFS.percent_good_within_LZA_threshold_qf;
    dqfflagstr11=num2str(dqfflag11,6);
    dqfflag12=100*DQFS.percent_degraded_due_to_LZA_threshold_exceeded_qf;
    dqfflagstr12=num2str(dqfflag12,6);
    br = PageBreak();
    add(chapter,br);
    DQFHdr={'Item' '% of Pixels'};
    DQFTable=cell(12,2);
    DQFTable{1,1}='Pct Overall Good Quality';
    DQFTable{1,2}=dqfflagstr1;
    DQFTable{2,1}='Pct Overall Degraded Quality';
    DQFTable{2,2}=dqfflagstr2;
    DQFTable{3,1}='Pct Good Quality L1b data';
    DQFTable{3,2}=dqfflagstr3;
    DQFTable{4,1}='Pct Degraded Quality L1b data';
    DQFTable{4,2}=dqfflagstr4;
    DQFTable{5,1}='Pct Good Quality Beta Ratio';
    DQFTable{5,2}=dqfflagstr5;
    DQFTable{6,1}='Pct Degraded Quality Beta Ratio';
    DQFTable{6,2}=dqfflagstr6;
    DQFTable{7,1}='Pct Ice Cloud Determ on Strong Rad Signal';
    DQFTable{7,2}=dqfflagstr7;
    DQFTable{8,1}='Pct Ice Cloud Determ on Weak Rad Signal';
    DQFTable{8,2}=dqfflagstr8;
    DQFTable{9,1}='Pct Good Quality Surface Emissivity';
    DQFTable{9,2}=dqfflagstr9;
    DQFTable{10,1}='Pct Degraded Quality Surface Emissivity';
    DQFTable{10,2}=dqfflagstr10;
    DQFTable{11,1}='Pct Good Within LZA Threshold';
    DQFTable{11,2}=dqfflagstr11;
    DQFTable{12,1}='Pct Degraded Due To LZA Threshold Exceedance';
    DQFTable{12,2}=dqfflagstr12;
    T4=[DQFHdr;DQFTable];
    tbl4=Table(T4);
    tbl4.Style = [tbl4.Style {Border('solid','black','3px')}];
    tbl4.HAlign='center';
    tbl4.TableEntriesHAlign = 'center';
    tbl4.ColSep = 'single';
    tbl4.RowSep = 'single';
    r = row(tbl4,1);
    r.Style = [r.Style {Color('red'),Bold(true)}];
    bt4 = BaseTable(tbl4);
    tabletitle = Text('DQF ABI L2 Cloud Top Phase');
    tabletitle.Bold = false;
    bt4.Title = tabletitle;
    bt4.TableWidth="7in";
    add(chapter,bt4);
% Add text for this table
    parastr21='The Data Quality Factor (DQF) table for the Cloud Top Phase Product is shown above.';
    parastr22='There are a total of 12 DQF factors-shown in the table above.';
    parastr23=strcat('The first DQF reported is to overall fraction of good estimates which is-',dqfflagstr1,'.');
    parastr24='Items 5 and 6 depend on the value a quantity called the beta ratio.';
    parastr25='The beta ratio is the ratio of the effective optical depth of the cloud in two channels most often 11 and 14 which are in the IR.';
    parastr26='The last two quality factors as is typical relation to the local zenith angle or LZA.';
    parastr29=strcat(parastr21,parastr22,parastr23,parastr24,parastr25,parastr26);
    p2 = Paragraph(parastr29);
    p2.Style = {OuterMargin("0pt", "0pt","10pt","10pt")};
   add(chapter,p2);
%   add(rpt,chapter);

end
% if(iframe==1)
%     fprintf(fid,'%s\n','------- Finished Plotting Cloud Top Phase Data ------');
% end
end


