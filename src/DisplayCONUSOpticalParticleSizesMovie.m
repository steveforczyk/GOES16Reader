function  DisplayCONUSOpticalParticleSizesMovie(titlestr,iframe)
% Display the Optical Cloud Particle Sizes on a map
% can add lightning to the map
% Note this routine is basically the same as DisplayCloudTopHeightsRev2.,
% Written By: Stephen Forczyk
% Created: Apr 29,2022
% Revised: -----
% Classification: Unclassified

global MetaDataS;
global PSDS;
global DQFS tS yS xS tBS goesImagerS yImgS yImgBS;
global xImgS xImgBS SatDataS GeoSpatialS ;
global PSDDayStatS PSDNightStatS OutlierPSDS;
global AlgoS GoesWaveBand;
global AlgoProdVersionContainerS ProcessParamS;
global QLZAS  QLZABS;
global TSZAS  TSZABS;
global DSZAS NSZAS;
global DRLZAS NRLZAS;
global DSZABS NSZABS;
global DRLZABS NRLZABS;
global RLZAS RSZAS;
global DayAlgoSZAS NightAlgoSZAS;
global PctDayPixelS PctNightPixelS PctTermPixelS;
global GRBErrorsS L0ErrorsS;
global westEdge eastEdge northEdge southEdge;
global GOESFileName;
global RptGenPresent iCreatePDFReport pdffilename rpt chapter;
global JpegCounter JpegFileList;
global DayMonthNonLeapYear DayMonthLeapYear CalendarFileName;
global DQFHdr DQFTable MiscHdr MiscTable;
global DQFCauses DQFNormed DQFLabels;
global RasterLats RasterLons FOVLOC MOVLOC agrid;
global GMTCTPFracS GMTCTPFracE GridMatFileName;
global FlashDataByFrame GroupPosByFrame itype iaddLightning;
global GroupMaxFlashDurations FileList GroupLongestFlashData ifirstCallGLMLightning;
global nc_filenames numncfiles ncpath;
global TempMovieName vTemp numMovieFrames nowMovieFrame;

global widd2 lend2;
global initialtimestr igrid ijpeg ilog imovie;
global vert1 hor1 widd lend;
global vert2 hor2 machine;
global chart_time;
global Fz1 Fz2;
global idirector mov izoom iwindow;
global matpath GOES16path;
global jpegpath fid;
global smhrpath excelpath ascpath;
global ipowerpoint PowerPointFile scaling stretching padding;
global ichartnum;
% additional paths needed for mapping
global matpath1 mappath;
global canadapath stateshapepath topopath;
global trajpath militarypath;
global figpath screencapturepath gridpath;
global shapepath2 countrypath countryshapepath usstateboundariespath;
global GOES16Band1path GOES16Band2path GOES16Band3path GOES16Band4path;
global GOES16Band5path GOES16Band6path GOES16Band7path GOES16Band8path
global GOES16Band9path GOES16Band10path GOES16Band11path GOES16Band12path;
global GOES16Band13path GOES16Band14path GOES16Band15path GOES16Band16path;
global GOES16CloudTopHeightpath GOES16BandPaths;
if((iCreatePDFReport==1) && (RptGenPresent==1))
    import mlreportgen.dom.*;
    import mlreportgen.report.*;
end
if(iframe==1)
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','-----Start plot routine Display CONUS OpticalParticleSizes-----');
end
eval(['cd ' gridpath(1:length(gridpath)-1)]);
% gridstr1=strcat('Loaded Map Grid From File=',GridMatFileName);
% fprintf(fid,'%s\n',gridstr1);
% load(GridMatFileName,'RasterLats','RasterLons');
% dispstr=strcat('Loaded Grid Raster File-',GridMatFileName);
% disp(dispstr)
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
    GridMatFileName='FDCloudParticleSizeBoundaries.mat';
    agrid=exist(GridMatFileName,'file'); 
    waitstr='Calculating Full Disk Geographic Raster';
elseif((iframe==1) && (itype==2))
    GridMatFileName='ConusCloudParticleSizeBoundaries.mat'; 
    agrid=exist(GridMatFileName,'file'); 
    waitstr='Calculating Conus Geographic Raster';
elseif((iframe==1) && (itype==3))
    GridMatFileName='Meso1CloudParticleSizeBoundaries.mat';
    agrid=exist(GridMatFileName,'file');
    waitstr='Calculating Meso1 Geographic Raster';
elseif((iframe==1) && (itype==4))
    GridMatFileName='Meso2CloudParticleSizeBoundaries.mat'; 
    agrid=exist(GridMatFileName,'file'); 
    waitstr='Calculating Meso2 Geographic Raster';
end
if((ishiftnow==1) && (itype==1))
    GridMatFileName='FDCloudParticleSizeBoundaries.mat';
    waitstr='Re-Calculating Full Disk Geographic Raster';  
    agrid=0;
elseif((ishiftnow==1) && (itype==2))
    GridMatFileName='ConusCloudParticleSizeBoundaries.mat';  
    waitstr='Re-Calculating Conus Geographic Raster'; 
    agrid=0;
elseif((ishiftnow==1) && (itype==3))
    GridMatFileName='Meso1CloudParticleSizeBoundaries.mat';
    waitstr='Re-Calculating Meso1 Geographic Raster'; 
    agrid=0;
elseif((ishiftnow==1) && (itype==4))
    GridMatFileName='Meso2CloudParticleSizeBoundaries.mat';
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
OptPS=PSDS.values;
[nrows,ncols]=size(OptPS);
ntotpixels=nrows*ncols;
OptPS1D=reshape(OptPS,nrows*ncols,1);
OptPS1DS=sort(OptPS1D);
a1=isnan(OptPS1DS);
numnan=sum(a1);
numvals2=(nrows*ncols)-(numnan+1);
OptPS1DSF=zeros(numvals2,1);
for k=1:numvals2
    OptPS1DSF(k,1)=OptPS1DS(k,1);
end
sgoodfrac=numvals2/(nrows*ncols);
% Get some data relating to the data collection
NumDayPixels=ntotpixels*PctDayPixelS.value;
NumNightPixels=ntotpixels*PctNightPixelS.value;
NumTwilightPixels=ntotpixels*PctTermPixelS.value;
NumOutlierDayPSDPixels=OutlierPSDS.values1;
MaxLZARetrievalAngle=QLZAS.value;
MaxSolarRetrievalAngle=DSZAS.value;
minDayPSDstr=num2str(PSDDayStatS.values1,6);
maxDayPSDstr=num2str(PSDDayStatS.values2,6);
meanDayPSDstr=num2str(PSDDayStatS.values3,6);
fracdaypixels=PctDayPixelS.value;
fracnightpixels=PctNightPixelS.value;
fractwilightpixels=PctTermPixelS.value;
NumDay=floor(NumDayPixels);
if(NumDay>0)
    fracoutlierdaypixels=OutlierPSDS.values1/NumDay;
    NumOutlierDay=floor(NumOutlierDayPSDPixels);
else
    fracoutlierdaypixels=0.0000;
    NumOutlierDay=0;
end
pixelstr1=strcat('Number of Day Pixels=',num2str(NumDay),'-fraction=',num2str(fracdaypixels,6));
fprintf(fid,'%s\n',pixelstr1);
NumNight=floor(NumNightPixels);
pixelstr2=strcat('Number of Night Pixels=',num2str(NumNight),'-fraction=',num2str(fracnightpixels,6));
fprintf(fid,'%s\n',pixelstr2);
NumTwilight=floor(NumTwilightPixels);
pixelstr3=strcat('Number of Twilight Pixels=',num2str(NumTwilight),'-fraction=',num2str(fractwilightpixels,6));
fprintf(fid,'%s\n',pixelstr3);
pixelstr4=strcat('Number of Outlier Day Pixels=',num2str(NumOutlierDay),'-fraction=',num2str(fracoutlierdaypixels,6));
fprintf(fid,'%s\n',pixelstr4);
PSDstr1=strcat('Min Day PSD=',minDayPSDstr,'-Mean Day PSD=',meanDayPSDstr,'-Max Day PSD=',maxDayPSDstr);
fprintf(fid,'%s\n',PSDstr1);
MiscTable=cell(9,3);
MiscTable{1,1}='Num Of Day Pixels';
MiscTable{1,2}=floor(NumDayPixels);
MiscTable{1,3}='-';
MiscTable{2,1}='Num Of Night Pixels';
MiscTable{2,2}=floor(NumNightPixels);
MiscTable{2,3}='-';
MiscTable{3,1}='Num Of Twilight Pixels';
MiscTable{3,2}=floor(NumTwilightPixels);
MiscTable{3,3}='-';
MiscTable{4,1}='Num Of Out Of Range Pixels';
MiscTable{4,2}=floor(NumOutlierDay);
MiscTable{4,3}='-';
MiscTable{5,1}='Max LZA Angle';
MiscTable{5,2}=floor(MaxLZARetrievalAngle);
MiscTable{5,3}='Angle In Degrees';
MiscTable{6,1}='Max Solar Retrieval Angle';
MiscTable{6,2}=floor(MaxSolarRetrievalAngle);
MiscTable{6,3}='Angle In Degrees';
MiscTable{7,1}='Min Daytime Particle Size';
MiscTable{7,2}=minDayPSDstr;
MiscTable{7,3}='microns';
MiscTable{8,1}='Mean Daytime Particle Size';
MiscTable{8,2}=meanDayPSDstr;
MiscTable{8,3}='microns';
MiscTable{9,1}='Max Daytime Particle Size';
MiscTable{9,2}=maxDayPSDstr;
MiscTable{9,3}='microns';
MiscHdr={'Item' 'Value' 'Units'};
% Fetch the map limits
westEdge=double(GeoSpatialS.geospatial_westbound_longitude);
eastEdge=double(GeoSpatialS.geospatial_eastbound_longitude);
southEdge=double(GeoSpatialS.geospatial_southbound_latitude);
northEdge=double(GeoSpatialS.geospatial_northbound_latitude);
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
% Get the scan start date and time data based on the file name
GMTCTPIS=60*MinuteS+SecondS;
GMTCTPFracS=HourS+GMTCTPIS/3600;
GMTCTPIE=60*MinuteE+SecondE;
GMTCTPFracE=HourE+GMTCTPIE/3600;
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
% axesm ('pcarree','Frame','on','Grid','on','MapLatLimit',[southEdge northEdge],...
%      'MapLonLimit',[westEdge eastEdge],'meridianlabel','on','parallellabel','on','plabellocation',10,'mlabellocation',20,...
%      'MLabelParallel','south');
set(gcf,'Position',[hor1 vert1 widd lend]);
set(gcf,'MenuBar','none');
cmapland=colormap('jet');
geoshow(RasterLats,RasterLons,OptPS,'DisplayType','texturemap');
demcmap(OptPS,256,[],cmapland);
hold on
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
title(titlestr);
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
pos=get(gca,'Position');
x0=pos(1,1);
y0=pos(1,2);
poscolor=get(hc,'Position');
cbarx0=poscolor(1,1)+poscolor(1,3)-.05;
cbary0=poscolor(1,2)+poscolor(1,4)+.02;
newaxesh=axes('Position',[0 0 1 1]);
set(newaxesh,'XLim',[0 1],'YLim',[0 1]);
tx1=.13;
ty1=y0-.03;
txtstr1=strcat('Start Scan-Y',num2str(YearS),'-D-',num2str(DayS),'-H-',...
    num2str(HourS),'-M-',num2str(MinuteS),'-S-',num2str(SecondS),'-----',MonthDayYearStr);
txt1=text(tx1,ty1,txtstr1,'FontWeight','bold','FontSize',12);
tx2=.80;
ty2=.76;
txtstr2='Particle Size-microns';
txt2=text(tx2,ty2,txtstr2,'FontWeight','bold','FontSize',12);
set(newaxesh,'Visible','Off');
% Save this chart
if(rm==1)
    figstr=strcat(titlestr,'.jpg');
    eval(['cd ' jpegpath(1:length(jpegpath)-1)]);
    actionstr='print';
    typestr='-djpeg';
    [cmdString]=MyStrcat2(actionstr,typestr,figstr);
    eval(cmdString);
    plotstr1=strcat('Saved Cloud Optical Particle Size Plots as-',figstr);
    disp(plotstr1);
    fprintf(fid,'%s\n',plotstr1);
end
pause(chart_time);
frame=getframe(gcf);
writeVideo(vTemp,frame);
close('all');
if((iCreatePDFReport==1) && (RptGenPresent==1))
    [ibreak]=strfind(GOESFileName,'_e');
    is=1;
    ie=ibreak(1)-1;
    GOESShortFileName=GOESFileName(is:ie);
    headingstr1=strcat('Cloud Optical Particle Size For-',GOESShortFileName);
    chapter = Chapter("Title",headingstr1);
    add(chapter,Section('Cloud Optical Particle Size Map'));
    imdata = imread(figstr);
    [nhigh,nwid,ndepth]=size(imdata);
    image = mlreportgen.report.FormalImage();
    image.Image = which(figstr);
    pdftxtstr=strcat('Cloud Optical Particle Size For File-',GOESShortFileName);
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
% Now add some text -start by decribing the COD Metric
    parastr1=strcat('The Data for this chart was taken from file-',GOESShortFileName,'.');
    parastr2='Plotted on the chart is the Cloud Optical Particle Size (PSD) which is  the diameter in microns.';
    parastr6=strcat(parastr1,parastr2);
    p1 = Paragraph(parastr6);
    p1.Style = {OuterMargin("0pt", "0pt","10pt","10pt")};
    add(chapter,p1);
% Description of cloud particle size
    parastr41=strcat('Cloud particle size can be estimated using data from the ABI instrument along with numerical weather models.');
    parastr42=strcat('The GOES weather satellites can use this statistic, along with Cloud Particle Size (CPS) to characterize the earth radiative energy budget.',...
        'In turn, the COD and CPS used together along with a number of climate models as powerful weather predictors.');
    parastr43='The ABI instrument can return this data at 2 km resolution at either 5 or 15 intervals in either a CONUS or full disk format.';
    parastr44='Limitations in the instrument and relatively week aerosol signals,the current ABI algorithm does not attempt to return values over bright areas.';
    parastr45='Examples of such areas include regions where sun glint is strong,and bright deserts to name a few.';
    parastr49=strcat(parastr41,parastr42,parastr43,parastr44,parastr45);
    p4 = Paragraph(parastr49);
    p4.Style = {OuterMargin("0pt", "0pt","10pt","10pt")};
    add(chapter,p4);
    dqfflag1=100*DQFS.percent_good_quality_qf;
    dqfflagstr1=num2str(dqfflag1,6);
    dqfflag2=100*DQFS.percent_day_algorithm_pixel_qf;
    dqfflagstr2=num2str(dqfflag2,6);
    dqfflag3=100*DQFS.percent_night_algorithm_pixel_qf;
    dqfflagstr3=num2str(dqfflag3,6);
    dqfflag4=100*DQFS.percent_degraded_quality_due_to_snow_or_sea_ice_qf;
    dqfflagstr4=num2str(dqfflag4,6);
    dqfflag5=100*DQFS.percent_degraded_quality_due_to_twilight_qf;
    dqfflagstr5=num2str(dqfflag5,6);
    dqfflag6=100*DQFS.percent_invalid_due_to_clear_conditions_qf;
    dqfflagstr6=num2str(dqfflag6,6);
    dqfflag7=100*DQFS.percent_invalid_due_LZA_threshold_exceeded_qf;
    dqfflagstr7=num2str(dqfflag7,6);
    dqfflag8=100*DQFS.percent_degraded_due_to_LZA_threshold_exceeded_qf;
    dqfflagstr8=num2str(dqfflag8,6);
    dqfflag9=100*DQFS.percent_invalid_due_to_not_geolocated_qf;
    dqfflagstr9=num2str(dqfflag9,6);
    dqfflag10=100*DQFS.percent_invalid_due_to_missing_or_bad_input_data_qf;
    dqfflagstr10=num2str(dqfflag10,6);
    dqfflag11=100*DQFS.percent_degraded_due_to_nonconvergence_qf;
    dqfflagstr11=num2str(dqfflag11,6);
    DQFCauses=zeros(8,1);
    DQFCauses(1,1)=dqfflag4;
    DQFCauses(2,1)=dqfflag5;
    DQFCauses(3,1)=dqfflag6;
    DQFCauses(4,1)=dqfflag7;
    DQFCauses(5,1)=dqfflag8;
    DQFCauses(6,1)=dqfflag9;
    DQFCauses(7,1)=dqfflag10;
    DQFCauses(8,1)=dqfflag11;
    sumcauses=sum(DQFCauses);
    DQFLabels=cell(1,8);
    DQFLabels{1,1}='Degraded due to snow/ice';
    DQFLabels{1,2}='Degraded due twlilight';
    DQFLabels{1,3}='Invalid due clear conditions';
    DQFLabels{1,4}='Invalid due LZA';
    DQFLabels{1,5}='Degraded due LZA';
    DQFLabels{1,6}='Invalid Due GeoLocation Errors';
    DQFLabels{1,7}='Invalid Due Bad Data';
    DQFLabels{1,8}='Invalid Due Algo NonConvergence';
    if(sumcauses>0)
        DQFNormed=DQFCauses*100/sumcauses;
    else
        DQFNormed=DQFCauses;
    end    
% Now build a DQF Table of key values
    DQFHdr = {'Item' '% Of Pixels'};
    DQFTable=cell(11,2);
    DQFTable{1,1}='Pct Good Quality Pixel';
    DQFTable{1,2}=dqfflagstr1;
    DQFTable{2,1}='Pct Day Algorithm Pixel';
    DQFTable{2,2}=dqfflagstr2;
    DQFTable{3,1}='Pct Night Algorithm Pixel';
    DQFTable{3,2}=dqfflagstr3;
    DQFTable{4,1}='Pct Degraded Pixels Snow/Ice';
    DQFTable{4,2}=dqfflagstr4;
    DQFTable{5,1}='Pct Degraded Pixels Twilight';
    DQFTable{5,2}=dqfflagstr5;
    DQFTable{6,1}='Pct Invalid Clear Conditions';
    DQFTable{6,2}=dqfflagstr6;
    DQFTable{7,1}='Pct Invalid LZA Threshold Exceeded';
    DQFTable{7,2}=dqfflagstr7;
    DQFTable{8,1}='Pct Degraded LZA Threshold Exceeded';
    DQFTable{8,2}=dqfflagstr8;
    DQFTable{9,1}='Pct Degraded Not Geolocated';
    DQFTable{9,2}=dqfflagstr9;
    DQFTable{10,1}='Pct Invalid Due To Bad Data';
    DQFTable{10,2}=dqfflagstr10;
    DQFTable{11,1}='Pct Degraded Due To NonConvergence';
    DQFTable{11,2}=dqfflagstr11;
    br = PageBreak();
    add(chapter,br);
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
    tabletitle = Text('DQF Table Cloud Optical Depth');
    tabletitle.Bold = false;
    bt4.Title = tabletitle;
    bt4.TableWidth="7in";
    add(chapter,bt4);
% Find which factor causes thre most invalid/degraded pixels
    DQFCauses=zeros(8,1);
    DQFCauses(1,1)=dqfflag4;
    DQFCauses(2,1)=dqfflag5;
    DQFCauses(3,1)=dqfflag6;
    DQFCauses(4,1)=dqfflag7;
    DQFCauses(5,1)=dqfflag8;
    DQFCauses(6,1)=dqfflag9;
    DQFCauses(7,1)=dqfflag10;
    DQFCauses(8,1)=dqfflag11;
    [DQFCausesSort,Index]=sort(DQFCauses,'descend');
    Index3=Index+3;
    biggest_reason_str=DQFTable{Index3(1,1),1};
    is=4;
    ie=length(biggest_reason_str);
    prime_cause_str=biggest_reason_str(is:ie);
    biggest_reason_value=DQFTable{Index3(1,1),2};
    dqfstr1=strcat('The biggest reason for rejected pixels is-',prime_cause_str,'-which results in-',biggest_reason_value,....
        '%-of all reporting pixels reporting invalid values.');
    fprintf(fid,'%s\n',dqfstr1);
% Add text for this table
    parastr11='The DQF table above in intended to inform the user regarding the factors that play into the quality of the data.';
    parastr12=strcat('In the first row of the table is the per centage of all pixels that returned good data-',dqfflagstr1,'%.');
    parastr13='The factors that reduce the number below 100% involve clear pixels,sea ice pixels,and pixels in night or with high zenith angles.';
    parastr14='Rows 4 through 11 provide a mode detailed breakout for invalid or degraded pixels.';
    parastr15=strcat('Inspection of the data in rows 4 through 11 reveal that the biggest cause for invalid pixels is-',prime_cause_str,'-which caused-',...
        biggest_reason_value,'-% of pixels to return inavlid values.','On the next page is a pie chart showing those DQF factors relating to invalid pixels.');
    parastr16=strcat('Another useful piece of data is to establish how much of the scene is daylit,for this frame of data it is-',dqfflagstr2,'% of pixels.');
    parastr19=strcat(parastr11,parastr12,parastr13,parastr14,parastr15,parastr16);
    p2 = Paragraph(parastr19);
    p2.Style = {OuterMargin("0pt", "0pt","10pt","10pt")};
    add(chapter,p2);
    CreateDQFPieChart(titlestr)
% Build a table with miscellaneous data related to the data collection
% limits
    T5=[MiscHdr;MiscTable];
    tbl5=Table(T5);
    tbl5.Style = [tbl5.Style {Border('solid','black','3px')}];
    tbl5.HAlign='center';
    tbl5.TableEntriesHAlign = 'center';
    tbl5.ColSep = 'single';
    tbl5.RowSep = 'single';
    r = row(tbl5,1);
    r.Style = [r.Style {Color('red'),Bold(true)}];
    bt5 = BaseTable(tbl5);
    tabletitle = Text('Miscellaneous Data Collection Parameters');
    tabletitle.Bold = false;
    bt5.Title = tabletitle;
    bt5.TableWidth="7in";
    add(chapter,bt5);
% add a descriptive paragraph for bt5
      parastr30='The table above provides additional details on this data collection.';
      parastr31='Inspection of the table shows that the first 3 rows detail how many pixels were viewing day,night or twilight regions of the earth.';
      parastr32='Row 4 shows how many pixels provide out of range estimates of cloud particle size-the range of allowable sizes extends from 0 to 100 microns.'; 
      parastr33='Rows 5 through 6 specify the limits for the Local Zenith and Solar Zenith Angles as they relate to measurement quality.';
      parastr34='The next 3 rows provide the minimum,mean and maximum values for the calculated cloud particle sizes.';
      parastr39=strcat(parastr30,parastr31,parastr32,parastr33,parastr34);
      p3 = Paragraph(parastr39);
      p3.Style = {OuterMargin("0pt", "0pt","10pt","10pt")};
      add(chapter,p3);
end
close('all');
fprintf(fid,'%s\n','-----Exit plot routine Display CONUS Optical Particle Size-----');
fprintf(fid,'\n');
end


