% This executive script was created to read a series of files to create
% a movie showing the ABI_L2_CMI Cloud Particel Size product over the course of up to 1 day
% at Full Disk,Conus and Meso levels
% Created: April 1,2022
% Written By: Stephen Forczyk
% Revised: -----

% Classification: Unclassified/Public Domain
%% Set Up Globals
global Datasets;
global BandDataS MetaDataS;
global CMIS DQFS tS yS xS tBS goesImagerS yImgS yImgBS;
global xImgS xImgBS SatDataS GeoSpatialS;
global ReflectanceS AlgoS ErrorS EarthSunS VersionContainerS;
global GoesWaveBand GoesWaveBand2 GoesWaveBandTable ESunS kappa0S PlanckS FocalPlaneS;
global GOESFileName MapFormFactor;
global GOES16BandPaths;
global isaveGOESLightningData isaveCMI;
global CountyBoundaryFile;
global CountyBoundaries StateFIPSFile;
global StateFIPSCodes NationalCountiesShp;
global USAStatesShapeFileList USAStatesFileName;
global GOESPlatformData itype;
global westEdgeL eastEdgeL northEdgeL southEdgeL ifovShift ifovShiftLimit;
global FOVLOC MOVLOC;
global MonthDayStr MonthDayYearStr DayMonthNonLeapYear DayMonthLeapYear;
global IDStr;

global UrbanAreasShapeFile NorthAmericaLakes FireSummaryFile;
global idebug isavefiles iDemFlag;
global iPrimeRoads iCountyRoads iCommonRoads iStateRecRoads iUSRoads iStateRoads;
global iLakes;
global NumProcFiles ProcFileList;
global RptGenPresent iCreatePDFReport pdffilename rpt chapter tocc ifirstCall;
global JpegCounter JpegFileList;
global ImageProcessPresent;
global iReport ifixedImagePaths iaddLightning LightningDataFile;
global nc_filenames numncfiles ncpath ;
global iSun iCityPlot USACityDataBase PlotCityDataBase CityPopLimit;
global TempMovieName vTemp numMovieFrames nowMovieFrame;
global HurricaneFile HurricaneName iHurricane HurricaneDataS MatFileName;
global GroupMaxFlashDurations LightFileList GroupLongestFlashData ifirstCallGLMLightning;
global FlashDataByFrame GroupPosByFrame;
global ioverrideGOESLimits westEdgeOveride eastEdgeOveride northEdgeOveride southEdgeOveride;
global pinc minc;
global PhaseHistory ResultsFile ResultsFileName;

global fid;
global widd2 lend2;
global initialtimestr igrid ijpeg ilog imovie;
global vert1 hor1 widd lend;
global vert2 hor2 machine;
global chart_time;
global Fz1 Fz2;
global idirector mov izoom iwindow;

global matpath GOES16path;
global jpegpath tiffpath;
global smhrpath excelpath ascpath;
global ipowerpoint PowerPointFile scaling stretching padding;
global ichartnum;
global ColorList RGBVals ColorList2 xkcdVals LandColors;
global orange bubblegum brown brightblue;
% additional paths needed for mapping
global matpath1 mappath matlabpath USshapefilepath dtedpath;
global moviepath;
global northamericalakespath logpath pdfpath govjpegpath;
global GOES16Band1path GOES16Band2path GOES16Band3path GOES16Band4path;
global GOES16Band5path GOES16Band6path GOES16Band7path GOES16Band8path
global GOES16Band9path GOES16Band10path GOES16Band11path GOES16Band12path;
global GOES16Band13path GOES16Band14path GOES16Band15path GOES16Band16path;
global GOES16CloudTopHeightpath GOES16CloudTopTemppath GOES16Lightningpath;
global GOES16ConusBand1path shapefilepath Countryshapepath figpath;
global gridpath countyshapepath nationalshapepath summarypath;
global DayMonthNonLeapYear DayMonthLeapYear CalendarFileName;
global baseimagepath selectedImageFolder;
global pwd;


clc;
%% Set Up Fixed Paths
% Set up some fixed paths for the data. Set the present working directory
% to where the this code is located
pret=which('ReadGoes16Datasets.m');
[islash]=strfind(pret,'\');
numslash=length(islash);
is=1;
ie=islash(2);
pwd=pret(is:ie);
matlabpath=strcat(pwd,'Matlab_Data\');
mappath='D:\Forczyk\Map_Data\Matlab_Maps\';
moviepath='D:\Goes16\Movies\';
dtedpath='F:\DTED\';
CalendarFileName='CalendarDays.mat';
shapefilepath='D:\Goes16\ShapeFiles\';
countyshapepath='D:\Forczyk\Map_Data\MAT_Files_Geographic\';
CountyBoundaryFile='CountyBoundingBoxes';
nationalshapepath='D:\Forczyk\Map_Data\NationalShapeFiles\';
NationalCountiesShp='cb_2018_us_county_500k.shp';
UrbanAreasShapeFile='cb_2018_us_ua10M_500k.shp';
USAStatesFileName='USAStatesShapeFileMapListRev4.mat';
USshapefilepath='D:\Forczyk\Map_Data\USStateShapeFiles\';
northamericalakespath='D:\Forczyk\Map_Data\Natural_Earth\Ten_Meter_Physical\';
NorthAmericaLakes='ne_10m_lakes_north_america.shp';
GOES16path='D:\Goes16\Imagery\July25_2020\Band01\';
GOES16Lightningpath='H:\GOES16\Matfiles\';
jpegpath='D:\Goes16\Imagery\Aug26_2020\Jpeg_Files\';
pdfpath='D:\Goes16\Imagery\Aug26_2020\PDF_Files\';
govjpegpath='D:\Goes16\Documents\Jpeg\';
figpath='D:\Goes16\Imagery\Oct_FigFiles\';
logpath='H:\GOES16\Imagery\Apr29_2020\Log_Files\';
summarypath='D:\Goes16\Summary_Files\';
tiffpath='D:\Forczyk\Map_Data\InterstateSigns\';
FireSummaryFile='SummaryFileData.mat';
baseimagepath='F:\GOES16\Imagery\';
Countryshapepath='D:\Forczyk\Map_Data\CountryShapefiles\';
gridpath='D:\Goes16\Grids\';
MatFileName='USACityData.mat';
LightningDataFile='Oct-22-2017-GroupLightningData-Frames-4320.mat';

%% Set up GOESPlatformData Cell that will hold data
GOESPlatformData=cell(1,1);
GOESPlatformData{1,1}='GOES Sat';
GOESPlatformData{1,2}='Frame';
GOESPlatformData{1,3}='itype';
GOESPlatformData{1,4}='westEdge';
GOESPlatformData{1,5}='eastEdge';
GOESPlatformData{1,6}='northEdge';
GOESPlatformData{1,7}='southEdge';
GOESPlatformData{1,8}='Frame Time UTC';
GOESPlatformData{1,9}='Lat change';
GOESPlatformData{1,10}='Lon change';
GOESPlatformData{1,11}='Recompute Grid';
GOESPlatformData{1,12}='GOES Lat';
GOESPlatformData{1,13}='GOES Lon';
GOESPlatformData{1,14}='File Name';
GOESPlatformData{1,15}='Spare 1';
%% Set up a global Cell to will output data on the phase
PhaseHistory=cell(1,1);
PhaseHistory{1,1}='frame #';
PhaseHistory{1,2}='File Name';
PhaseHistory{1,3}='Frame Start Time';
PhaseHistory{1,4}='Frame End Time';
PhaseHistory{1,5}='NaN Fraction';
PhaseHistory{1,6}='Clear Sky Fraction';
PhaseHistory{1,7}='Liquid Water Fraction';
PhaseHistory{1,8}='Super Cooled Water Fraction';
PhaseHistory{1,9}='Mixed Fraction';
PhaseHistory{1,10}='Ice Fraction';
PhaseHistory{1,11}='Unknown Fraction';

%%
%% Set Flags
% Set some flags to control program execution and initialize a few counters
iCreatePDFReport=0;
JpegCounter=0;
isavefiles=0;
idebug=0;
ifixedImagePaths=1;
iSun=1;
iHurricane=0;
iCityPlot=0;
CityPopLimit=500000;
iaddLightning=0;
ifirstCall=0;
%% Set up a global cell to hold the location of the field of view data.
% when this changes by a sufficient amount these value will trigger a
% recompute of the Raster Grid data. This will likely occur only in meso
% data
ifovShift=0;
ifovShiftLimit=1.0;
FOVLOC=cell(1,1);
FOVLOC{1,1}='Frame #';
FOVLOC{1,2}='GMT Start Time';
FOVLOC{1,3}='Existing West Edge';
FOVLOC{1,4}='Existing East Edge';
FOVLOC{1,5}='Existing North Edge';
FOVLOC{1,6}='Existing South Edge';
FOVLOC{1,7}='New West Edge';
FOVLOC{1,8}='New East Edge';
FOVLOC{1,9}='New North Edge';
FOVLOC{1,10}='New South Edge';
FOVLOC{1,11}='4 Corner Movement Deg';
FOVLOC{1,12}='Max Allowed Movement Deg';
FOVLOC{1,13}='Shift Required';
FOVLOC{1,14}='Shift #';
MOVLOC=cell(1,1);
MOVLOC=FOVLOC;
%% Set up some starting values
% The next 4 values are just dummy starting values
westEdgeL=-300;
eastEdgeL=-280;
northEdgeL=87;
southEdgeL=84;
%% Set Up Log File
% Start writing a log file and Also look at the current stored image paths
% file
%eval(['cd ' logpath(1:length(logpath)-1)]);
startruntime=deblank(datestr(now));
startrunstr=strcat('Start Run GOES 16 Run at-',startruntime);
logfilename=startruntime;
logfiledbl=double(logfilename);
% find the blank space in the date and replace it with a '-' to make a
% legalfilename
[iblank]=find(logfiledbl==32);
numblank=length(iblank);
for n=1:numblank
    is=iblank(n);
    ie=is;
    logfilename(is:ie)='-';
end
[icolon]=strfind(logfilename,':');
numcolon=length(icolon);
for n=1:numcolon
    is=icolon(n);
    ie=is;
    logfilename(is:ie)='-';
end
datetimestr=logfilename;
ImageFileFoldersName=strcat('ImageFolderFiles-',datetimestr,'.mat');
[idash]=strfind(ImageFileFoldersName,'-');
numdash=length(idash);
[iper]=strfind(ImageFileFoldersName,'.');
numper=length(iper);
is=idash(1)+1;
ie=iper(1)-1;
PotentialDatestr=ImageFileFoldersName(is:ie);
idash=strfind(PotentialDatestr,'-');
is=idash(3);
ie=is;
PotentialDatestr(is:ie)=' ';
is=idash(4);
ie=is;
PotentialDatestr(is:ie)=':';
is=idash(5);
ie=is;
PotentialDatestr(is:ie)=':';
pd=datenum(PotentialDatestr,'dd-mmm-yyyy HH:MM:SS');
dirlis=dir(matlabpath);
% Find out how many files start with 'ImageFolderFiles'
numdirfiles=length(dirlis)-2;
maxage=0;
minage=100;
ibest=0;
if(numdirfiles<1)
    iwrite=0;
else
    pattern='ImageFolderFiles';
    FoundFiles=cell(1,4);
    ihit=0;
    ihitbest=1;
    for kk=3:numdirfiles+2
        nowFile=dirlis(kk).name;
        nowDateNum=dirlis(kk).datenum;
        TF = startsWith(nowFile,pattern);
        if(TF==1)
            ihit=ihit+1;
            FoundFiles{ihit,1}=nowFile;
            FoundFiles{ihit,2}=nowDateNum;
            FoundFiles{ihit,3}=pd-nowDateNum;
            FoundFiles{ihit,4}=0;% Keep flag
            nowage=pd-nowDateNum;
            if(nowage>maxage)
                maxage=nowage;
            end
            if(nowage<minage)
                minage=nowage;
                ibest=kk;
                ihitbest=ihit;
            end
        end
    end
    FoundFiles{ihitbest,4}=1;
    if(maxage>1)
        iwrite=1;
    else
        iwrite=0;
    end
    if(ihit>1)
    eval(['cd ' matlabpath(1:length(matlabpath)-1)]);
        for jj=1:ihit
            nowFile=strcat(matlabpath,char(FoundFiles{jj,1}));
            nowFlag=FoundFiles{jj,4};
            if(nowFlag==0)
                delete(nowFile);
            end
        end
    end
end

%% Set file source if default Image Path or a New Image File path will be set up
if((ifixedImagePaths>2) && (iwrite==1))
    ipathstr1='User has selected the option to set new image data filepaths';
    fprintf(fid,'%s\n',ipathstr1);
    selectedImageFolder=[];
    pret=which('ReadGoes16Datasets.m');
    [islash]=strfind(pret,'\');
    numslash=length(islash);
    is=1;
    ie=islash(numslash);
    pwd=pret(is:ie);
    [inFolderList] = GetSubDirsFirstLevelOnly(baseimagepath);
    [outFolderList] = MakeDateFolderList(baseimagepath,inFolderList);
% Save this list for future use
    eval(['cd ' matlabpath(1:length(matlabpath)-1)]);
    datetimestr2=datetime('now');
    dateID=date;
    actionstr='save';
    varstr='outFolderList dateID';
    qualstr='-v7.3';
    [cmdString]=MyStrcatV73(actionstr,ImageFileFoldersName,varstr,qualstr);
    eval(cmdString)
    dispstr=strcat('Wrote Matlab File-',ImageFileFoldersname);
    disp(dispstr);
    fprintf(fid,'%s\n',dispstr);
    [indx,tf] = listdlg('PromptString',{'Select folder to read'},...
        'SelectionMode','single','ListString',outFolderList,'ListSize',[360,300]);
    a1=isempty(indx);
    if(a1==1)
        ifixedImagePaths=1;
        ipathstr2='User did not select a new filepath-go to default paths';
        fprintf(fid,'%s\n',ipathstr2);
    else
        selectedImageFolder=outFolderList{indx,1};
        ipathstr3=strcat('User selected folder-',selectedImageFolder,'-as the source of image data');
        fprintf(fid,'%s\n',ipathstr3);
    end
elseif((ifixedImagePaths>2) && (iwrite==0))
    eval(['cd ' matlabpath(1:length(matlabpath)-1)]);
    loadfile=dirlis(ibest).name;
    ipathstr2='User has selected the option to select new image data filepaths but no need to recompute available folderlist';
    fprintf(fid,'%s\n',ipathstr2);
    load(loadfile);
    [indx,tf] = listdlg('PromptString',{'Select prev folder to read'},...
        'SelectionMode','single','ListString',outFolderList,'ListSize',[360,300]);
    a1=isempty(indx);
    if(a1==1)
        ifixedImagePaths=1;
        ipathstr2='User did not select a new filepath-go to default paths';
        fprintf(fid,'%s\n',ipathstr2);
    else
        selectedImageFolder=outFolderList{indx,1};
        ipathstr3=strcat('User selected folder-',selectedImageFolder,'-as the source of image data');
        fprintf(fid,'%s\n',ipathstr3);
    end
    
end
%% User chose to go with fixed "default file path" set these up
if(ifixedImagePaths==1)
    GOES16BandPaths=cell(77,1);
    GOES16BandPaths{1,1}='H:\GOES16\Imagery\Apr29_2020\ABI_L2_Cloud_Particle_Size\Full_Disk\';
    GOES16BandPaths{2,1}='H:\GOES16\Imagery\Apr29_2020\ABI_L2_Cloud_Particle_Size\Conus\';
    GOES16BandPaths{3,1}='H:\GOES16\Imagery\Apr29_2020\ABI_L2_Cloud_Particle_Size\Meso1\';
    GOES16BandPaths{4,1}='H:\GOES16\Imagery\Apr29_2020\ABI_L2_Cloud_Particle_Size\Meso2\';
    jpegpath='H:\GOES16\Imagery\Apr29_2020\Jpeg_Files\';
    pdfpath='H:\GOES16\Imagery\Apr29_2020\PDF_Files\';
    logpath='H:\GOES16\Imagery\Apr29_2020\Log_Files\';
    moviepath='H:\GOES16\Imagery\Apr29_2020\Movies\';
    eval(['cd ' logpath(1:length(logpath)-1)]);
    logfilename=strcat('CloudParticleSizeLogFile-',logfilename,'.txt');
    pdffilename=strcat('GOES16-',datetimestr);
    fid=fopen(logfilename,'w');
    dispstr=strcat('Opened Log file-',logfilename,'-for writing');
    disp(dispstr);
    fprintf(fid,'%50s\n',startrunstr);
    ipathstr1='User has selected the option to go with default image file paths used in distribution';
    fprintf(fid,'%s\n',ipathstr1);
    LightningDataFile='Apr-29-2020-GroupLightningData-Frames-4320.mat';
    ioverrideGOESLimits=0;
    westEdgeOveride = -100;
    eastEdgeOveride = -80;
    northEdgeOveride = 40;
    southEdgeOveride = 20;
    pinc=2;
    minc=2;
    ResultsFile='CloudTopData-Apr29-2020';
    IDStr='Apr29-2020';
%% user chose to run with data repository #2
elseif(ifixedImagePaths==2)% This choice is to be used in testing files that will be part of the distribution
    GOES16BandPaths=cell(1,1);
    GOES16BandPaths{1,1}='H:\GOES16\Imagery\Oct22_2017\ABI_L2_Cloud_Particle_Size\Full_Disk\';
    GOES16BandPaths{2,1}='H:\GOES16\Imagery\Oct22_2017\ABI_L2_Cloud_Particle_Size\Conus\';
    GOES16BandPaths{3,1}='H:\GOES16\Imagery\Oct22_2017\ABI_L2_Cloud_Particle_Size\Meso1\';
    GOES16BandPaths{4,1}='H:\GOES16\Imagery\Oct22_2017\ABI_L2_Cloud_Particle_Size\Meso1\';
    jpegpath='H:\GOES16\Imagery\Oct22_2017\Jpeg_Files\';
    pdfpath='H:\GOES16\Imagery\Oct22_2017\PDF_Files\';
    logpath='H:\GOES16\Imagery\Oct22_2017\Log_Files\';
    moviepath='H:\GOES16\Imagery\Oct22_2017\Movies\';
    eval(['cd ' logpath(1:length(logpath)-1)]);
    logfilename=strcat('CloudParticleSizeLogFile-',logfilename,'.txt');
    pdffilename=strcat('GOES16-',datetimestr);
    fid=fopen(logfilename,'w');
    dispstr=strcat('Opened Log file-',logfilename,'-for writing');
    disp(dispstr);
    fprintf(fid,'%50s\n',startrunstr);
    ipathstr1='User has selected the option to go with default image file paths used in distribution';
    fprintf(fid,'%s\n',ipathstr1);
    LightningDataFile='Oct-22-2017-GroupLightningData-Frames-4320.mat';
    ioverrideGOESLimits=0;
    westEdgeOveride = -100;
    eastEdgeOveride = -80;
    northEdgeOveride = 40;
    southEdgeOveride = 20;
    pinc=2;
    minc=2;
    westEdgeOveride = -100;
    eastEdgeOveride = -94;
    northEdgeOveride = 38;
    southEdgeOveride = 34;
    pinc=0.5;
    minc=0.5;
    ResultsFile='CloudTopData-Oct22-2017';
    IDStr='Oct22-2017';
else
%% User chose to add a new Image File data (different day) make these paths
% Now set the image paths to the selected folder adding any missing
% subfolders as required
    [output1,output2] = SetImageFolders();
end
MapFormFactor=[];
%% Set up Dateset Names
Datasets=cell(1,1);
Datasets{1,1}='ABI-L2-Cloud Particle Size Full Disk';
Datasets{2,1}='ABI-L2-Cloud Particle Size Conus';
Datasets{3,1}='ABI-L2-Cloud Particle Size Meso 1';
Datasets{4,1}='ABI-L2-Cloud Particle Size Meso 2';
Datasets{5,1}='Stop Run';
iReport=zeros(78,1);

%% Set up some initial data
ProcFileList=cell(1,3);
ProcFileList{1,1}='File Type';
ProcFileList{1,2}='File Name';
ProcFileList{1,3}='Source Directory';
JpegFileList=cell(1,2);
JpegFileList{1,1}='Source Directory';
JpegFileList{1,2}='JpegFileName';
SetUpGoesWaveBandData();
%% Check to see if the Matlab Report Generator and Image Toolbox is present
toolbox='MATLAB Report Generator';
[RptGenPresent] = ToolboxChecker(toolbox);
toolbox='Image Processing Toolbox';
[ImageProcessPresent]=ToolboxChecker(toolbox);
%% Start a PDF report if requested by user and user has the Matlab Report
% Generator Package Installed
if((iCreatePDFReport==1) && (RptGenPresent==1))
    import mlreportgen.dom.*;
    import mlreportgen.report.*;
    import mlreportgen.utils.*
    CreatePDFReportHeaderRev7
else
    fprintf(fid,'%s\n','No PDF report will be created for this run');
end
%% Read Some needed data files related to calendar data and selected shape file names
eval(['cd ' matlabpath(1:length(matlabpath)-1)]);
load(CalendarFileName);
load('StandardAtmosphere.mat');
% Load in the CountyBoundaryFiles
eval(['cd ' countyshapepath(1:length(countyshapepath)-1)]);
load('CountyBoundingBoxes.mat');
% Load in the list of USAStateShapeFiles
load(USAStatesFileName);
% Set up some colors that will be used in later plots
SetUpExtraColors()
%% Call some routines that will create nice plot window sizes and locations
% Establish selected run parameters
imachine=2;
if(imachine==1)
    widd=720;
    lend=580;
    widd2=1000;
    lend2=700;
elseif(imachine==2)
    widd=1080;
    lend=812;
    widd2=1000;
    lend2=700;
elseif(imachine==3)
    widd=1296;
    lend=974;
    widd2=1200;
    lend2=840;
end
% Set a specific color order
set(0,'DefaultAxesColorOrder',[1 0 0;
    1 1 0;0 1 0;0 0 1;0.75 0.50 0.25;
    0.5 0.75 0.25; 0.25 1 0.25;0 .50 .75]);
% Set up some defaults for a PowerPoint presentationwhos
scaling='true';
stretching='false';
padding=[75 75 75 75];
igrid=1;
% Set up paramters for graphs that will center them on the screen
[hor1,vert1,Fz1,Fz2,machine]=SetScreenCoordinates(widd,lend);
[hor2,vert2,Fz1,Fz2,machine]=SetScreenCoordinates(widd2,lend2);
chart_time=5;
idirector=1;
initialtimestr=datestr(now);
igo=1;
NumProcFiles=0;
%% Read in Lightning Data from a previous run of the script
% CreateLightningMovie if it is desired to add lightning data to the CMI
% moisture plot
if(iaddLightning>0)
    eval(['cd ' GOES16Lightningpath(1:length(GOES16Lightningpath)-1)]);
    result = isfile(LightningDataFile);
    if(result==1)
        load(LightningDataFile,'FileList','FlashDataByFrame','GroupLongestFlashData','GroupMaxFlashDurations','GroupPosByFrame');
        LightFileList=FileList;
        dispstr=strcat('Successfully loaded Lightning Data from File-',LightningDataFile);
        disp(dispstr);
        lightstr=strcat('Loaded file-',LightningDataFile,'-will be used to add lightning data to movie');
        fprintf(fid,'%s\n',lightstr);
    else
        igo=0;
        disp('Sorry had to terminate run as lightning file could not be found')
    end
end

iframe=0; % Movie frame number
%% This is the main executive loop of the routine
while igo>0 % This setup up a loop to processing various file until user decides to STOP RUN
    [indx,tf] = listdlg('PromptString',{'Select type of data to read'},...
    'SelectionMode','single','ListString',Datasets,'ListSize',[360,300]);
    a1=isempty(indx);
    if(a1==1)
        igo=0;
    else
        igo=1;
    end
    if(indx==5)
        igo=0;
    end
    if(igo==0)
        break
    end
    if(indx==1)
        MapFormFactor='Full-Disk';
        itype=1;
    elseif(indx==2)
        MapFormFactor='Conus';
        itype=2;
    elseif(indx==3)
        MapFormFactor='Meso1';
        itype=3;
    elseif(indx==4)
        MapFormFactor='Meso2';
        itype=4;
    end
GOES16path=GOES16BandPaths{indx,1};
% Go to the expected path
    eval(['cd ' GOES16path(1:length(GOES16path)-1)]);
    if(indx<5)
        [nc_filenames,ncpath]=uigetfile('*nc','Select First and Last file GOES Data File','MultiSelect','on');
        numncfiles=length(nc_filenames);
        nc_filenames=nc_filenames';
        iframe=iframe+1;
        numMovieFrames=numncfiles;
        GOESFileName1=nc_filenames{iframe,1};
        GOESFileName2=nc_filenames{numncfiles,1};
        GOESFileName=GOESFileName1;
        numfiles=length(nc_filenames);
        MakeOpticalSizeMovieNameFromGOESDataFile(numfiles)
        eval(['cd ' moviepath(1:length(moviepath)-1)]);
        vTemp = VideoWriter(TempMovieName,'MPEG-4');
        vTemp.Quality=100;
        if(indx<3)
            vTemp.FrameRate=5;
        elseif(indx>2)
            vTemp.FrameRate=20;
        end
        open(vTemp);
        dispstr=strcat('Opened Video Object to create movie=',TempMovieName);
        disp(dispstr)
        for iframe=1:numncfiles
            ReadCloudParticleSizeMov(iframe)
            dispstr=strcat('Captured Frame-',num2str(iframe),'-of-',num2str(numncfiles),...
                '-Frames');
            disp(dispstr)
            if(iframe>1)
                idebug=0;
            end
        end
        close(vTemp);
        dispstr=strcat('Finished creating movie=',TempMovieName,'-which had-',num2str(numMovieFrames),'-Frames');
        disp(dispstr)
        igo=0;
        if((iCreatePDFReport==1) && (RptGenPresent==1))
            add(rpt,chapter);
        end

    end
end

%% Save the run results
ResultsFileName=strcat(ResultsFile,'-',MapFormFactor,'-FR-',num2str(numfiles),'.mat');
% eval(['cd ' savepath(1:length(savepath)-1)]);
% actionstr='save';
% varstr='PhaseHistory FOVLOC MOVLOC';
% qualstr='-v7.3';
% [cmdString]=MyStrcatV73(actionstr,ResultsFileName,varstr,qualstr);
% eval(cmdString)
% dispstr=strcat('Wrote Matlab File-',ResultsFileName);
% disp(dispstr);
% fprintf(fid,'\n');
% savestr1=strcat('User saved Cloud Top Phase Data to file-',ResultsFileName);
% fprintf(fid,'%s\n',savestr1);
% Plot a Summary of the cloud type fractions over the run
if(itype==1)
    titlestr='CloudTopRun-Full-Disk-Apr29-2020';
elseif(itype==2)
    titlestr='CloudTopRun-Conus-Apr29-2020';
elseif(itype==3)
    titlestr='CloudTopRun-Meso1-Apr29-2020';
elseif(itype==4)
    titlestr='CloudTopRun-Meso2-Apr29-2020';
end
%PlotCloudTopPhaseFractions(titlestr)
% Now plot the FOV Shift Changes
if(itype==1)
    titlestr='CloudTopRunFOVChanges-Full-Disk-Apr29-2020';
elseif(itype==2)
    titlestr='CloudTopRunFOVChanges-Conus-Apr29-2020';
elseif(itype==3)
    titlestr='CloudTopRunFOVChanges-Meso1-Apr29-2020';
elseif(itype==4)
    titlestr='CloudTopRunFOVChanges-Meso2-Apr29-2020';
end
PlotFOVChangeTimes(titlestr)
%% Run closeout
endruntime=deblank(datestr(now));
endrunstr=strcat('Finished GOES 16 Run at-',endruntime);
fprintf(fid,'%s\n',endrunstr);
fclose(fid);
% Close a PDF report if one is created
a1=exist('rpt','var');
if((iCreatePDFReport==1) && (RptGenPresent==1))
    close(rpt);
    rptview(rpt)
    dispstr=strcat('Closed PDF Report-',pdffilename);
    disp(dispstr)
else
    disp('No pdf report generated by this run');
end

ab=1;
