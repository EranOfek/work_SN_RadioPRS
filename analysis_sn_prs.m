%%

% downloaded on 6-Mar-2026
[TNS]=VO.TNS.downloadAll;
size(TNS,1)
% 190,894

ColNames = {'objid','name','ra','declination','type','redshift','discoverymag','discoverydate'};

TNS = TNS(:,ColNames);
TNS.Origin = zeros(size(TNS,1),1);


% merge with TNS_old
load('TNS_old.mat'); % Load the old TNS data

TNS_old = TNS_old(:,{'ID','Name','RA','DEC','Obj_Type','Redshift','DiscoveryMag_Flux','DiscoveryDate_UT_'});
TNS_old.Properties.VariableNames = ColNames;
TNS_old.Origin = ones(size(TNS_old,1),1);

TNS = [TNS;TNS_old];




% 7-3-2026 (17:18)
% remove FRBs
FF= TNS.type~="FRB";
TNS = TNS(FF,:);
size(TNS,1)
% left with 185,845

%% add redshifts
R=imUtil.cat.match2Galaxies(TNS.ra,TNS.declination);

TNS = tools.table.addStructFields2TableCols(TNS, R);
    
Z = tools.array.selectFirstNotNaN(TNS.redshift, TNS.Z_PGC, TNS.Z_GLADE);

TNS.Z = Z;

%% save table
save -v7.3 TNS.mat TNS


  
%%

load TNS.mat;



% write table of TNS with redshifts
AllWithZ  = TNS(:,{'objid','name','type','ra','declination','redshift','Z'});
tools.table.sprintf_table(AllWithZ(1:5,:),'Format',{'%7d','%-7s','%-9s','%10.6f','%10.6f','%6.4f','%6.4f'},'IsLatex',true)
TableStr=tools.table.sprintf_table(AllWithZ,'Format',{'%7d','%-7s','%-9s','%10.6f','%10.6f','%6.4f','%6.4f'},'IsLatex',false);
FID = fopen('Table_TNS_WithZ.txt','w');
fprintf(FID,'%s',TableStr);
fclose(FID);

%%

FlagZ=~isnan(TNS.Z); % <0.1; %06;    
TNS=TNS(FlagZ,:);

% 17,397 (before adding Z) with z<0.1
% 53,266 (after adding Z) with z<0.1
% 98,676 non-nan Z

%%
RAD = 180./pi;

SearchRadius = 10;  % [arcsec]

Nsn = size(TNS,1);
tic;
for Isn=1:1:Nsn
	  [Isn, Nsn]
    [CatRadio, ColRadio, ~, Dist] = catsHTM.cone_search('VLASSep1', TNS.ra(Isn)./RAD, TNS.declination(Isn)./RAD, SearchRadius);
     
    if Isn==1
        Ncol = numel(ColRadio);
        MatchedRadio = nan(Nsn, Ncol);
        MatchedDist  = [nan(Nsn, 1), zeros(Nsn,1)];
    end
    if ~isempty(CatRadio)
        Dist = Dist.* RAD.*3600;
        [MinDist, Imin] = min(Dist);
        CatRadio = CatRadio(Imin,:);

        MatchedRadio(Isn,:) = CatRadio;
        MatchedDist(Isn,:)  = [MinDist, numel(Dist)];
    end
end
toc

save -v7.3 Matched.mat MatchedRadio MatchedDist

sum(MatchedDist(:,1)<1)
% 1342 unique matches

    
%%
% select radio point sources
Fps = ~isnan(MatchedRadio(:,5)) & MatchedRadio(:,5)< (MatchedRadio(:,7) + 3.* sqrt(MatchedRadio(:,6).^2 + MatchedRadio(:,8).^2));
sum(Fps)
% 2236 selected

VLASS = cats.radio.VLASS1;
FPS=~isnan(VLASS.Table.Ftot) & VLASS.Table.Ftot<(VLASS.Table.Fpeak + 3.*sqrt(VLASS.Table.ErrFtot.^2 + VLASS.Table.ErrFpeak));
sum(FPS)./size(VLASS.Catalog,1)
% 81% are point sources

% select clear associations (<1"):
Fd=MatchedDist(:,1)<1; 
% 1343 unique matches

% redshift
Fz = TNS.Z<0.1;
% 53,266 matches

% discovered prior to 2017 (VLASS epoch 1):
Ftime = TNS.discoverydate.Year<2017;
% 3973 matches

Fprs = Fps & Fd & Fz & Ftime;
sum(Fprs)
% 10 matches


% abs mag:
TNS.disc_abs_mag = TNS.discoverymag - (5.*log10(astro.cosmo.lum_dist(TNS.Z))-5);

save -v7.3 TNS_Selected.mat TNS

%%

Selected = TNS(Fprs,{'name','type','ra','declination','Z','discoverymag','internal_names','disc_abs_mag'});
Selected = TNS(Fprs,{'name','type','ra','declination','Z','discoverymag','disc_abs_mag'});

% match to FIRST
FIRST=cats.radio.FIRST;
Ncand = size(Selected,1);
Selected.FIRST = nan(Ncand,1);
Selected.DistFIRST = nan(Ncand,1);
NVSS = cats.radio.NVSS;
NVSS_Dist = nan(Ncand,1);
Selected.NVSS = nan(Ncand,1);
for Icand=1:1:Ncand
    D = celestial.coo.sphere_dist_fast(Selected.ra(Icand)./RAD, Selected.declination(Icand)./RAD,FIRST.Cat(:,1),FIRST.Cat(:,2));
    [DistFIRST, IndFIRST] = min(D);
    Selected.DistFIRST(Icand) = DistFIRST.*RAD.*3600;
    if (DistFIRST.*RAD.*3600)<1.5
        Selected.FIRST(Icand) = FIRST.Cat(IndFIRST,4);  % Fpeak [mJy]

    end 
    D = celestial.coo.sphere_dist_fast(Selected.ra(Icand)./RAD, Selected.declination(Icand)./RAD,NVSS.Cat(:,1),NVSS.Cat(:,2));
    [NVSS_Dist(Icand), Invss] = min(D);
    NVSS_Dist(Icand) = NVSS_Dist(Icand).*RAD.*3600;
    if NVSS_Dist(Icand)<10
        Selected.NVSS(Icand) = NVSS.Cat(Invss,5);
    end
end
Selected.DistFIRST = Selected.DistFIRST<1000;  % in FIRST footprint

SelectedRadio = array2table(MatchedRadio(Fprs,:),'VariableNames',VLASS.ColNames);
SelectedRadio = SelectedRadio(:,{'Fpeak', 'ErrFpeak'});
SelectedRadio.nuLnu = 1e-26.*3e9.*SelectedRadio.Fpeak.*4.*pi.*(astro.cosmo.lum_dist(Selected.Z).*constant.pc).^2;
tools.table.sprintf_table([Selected, SelectedRadio],'Format',{'%-7s','%-9s','%10.6f','%10.6f','%6.4f','%4.1f','%6.1f','%5.2f','%1d','%5.2f','%5.2f','%5.2f','%4.1e'},'IsLatex',true)
[Selected, SelectedRadio]


% select transients with no NVSS/FIRST detection:
Ind_Fprs = find(Fprs);
Ind_Fprs = Ind_Fprs(~(~isnan(Selected.FIRST) | ~isnan(Selected.NVSS)))

%% plots and rate

% fraction of transients detected in 2016
sum(TNS.discoverydate.Year==2016)./sum(TNS.discoverydate.Year<2017)

% z distribution of all TNS prior to 2017 and have z
F=TNS.discoverydate.Year<2017;  
%TNS=TNS(F,:);

Zed=(0:0.01:1);
[Nz]=histcounts(TNS.Z(F),Zed);       
Xz=(Zed(1:end-1)+Zed(2:end))./2;
bar(Xz,Nz) 


% compare the histogram of redshift with the transients with z distribution
%[Nprsz]=histcounts(Selected.Z, Zed);       

%bar(Xz, Nprsz./Nz)




% histogram of discvery year
YearEd = (2000:1:2017).';
Xyear  = (YearEd(1:end-1)+YearEd(2:end)).*0.5;
[Nyear]=histcounts(TNS.discoverydate.Year,YearEd)
bar(Xyear, Nyear)

%%
F = TNS.discoverydate.Year<2017;  
[~,Volume] = astro.cosmo.comoving_volume(TNS.Z(Ind_Fprs));
MeanVolume = mean(Volume);
[MeanVolumeZ]=astro.cosmo.inv_comoving_volume(MeanVolume);
fprintf('z of mean volume: %f\n',MeanVolumeZ);
DistMeanVol   = astro.cosmo.lum_dist(MeanVolumeZ);
fprintf('Dist of mean volume: %e\n',DistMeanVol);
fprintf('nuLnu of mean volume distance: %e\n',3e9.*1.*1e-26.*4.*pi.*(DistMeanVol.*constant.pc).^2);

% rate compared to SN rate:
fprintf('Number of transients with mean volume z: %d\n',sum(F & TNS.Z<MeanVolumeZ));
fprintf('Rate in units of SN rate: %f\n',numel(Ind_Fprs).*0.5./sum(F & TNS.Z<MeanVolumeZ));

%% fit rate and lum. fun.
% 3D fit: lum function cutoff, slope, and rate

F=TNS.discoverydate.Year<2017;  

Zed=(0:0.001:1);
[Nz]=histcounts(TNS.Z(F),Zed);   
Xz=(Zed(1:end-1)+Zed(2:end))./2;
bar(Xz,Nz) 

VecLpeak = logspace(35,40,100);
VecSlope = (-1:0.5:2).';
N_Lpeak  = numel(VecLpeak);
N_Slope  = numel(VecSlope);

VecL = logspace(35,40,100).';
for Ilp=1:1:N_Lpeak
    for Isl=1:1:N_Slope


        Xpeak = VecLpeak(Ilp);
        Slope = VecSlope(Isl);

        Par    = [Xpeak, Npeak, Slope];
        LumFun = tools.math.fun.funPowerLawExpCutoff(Par, VecL);
    end
end



% for each SN in list TNS - calc its Dist:
DistTNS = astro.cosmo.lum_dist(TNS.Z).*constant.pc;  % [cm]
Dist = astro.cosmo.lum_dist(Xz).*constant.pc;

