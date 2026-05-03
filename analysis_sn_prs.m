%%

% downloaded on 3-May-2026  07:55 UTC
[TNS]=VO.TNS.downloadAll('AddOld',true);
size(TNS,1)
% 200,855

%% save table
save -v7.3 TNS_All.mat TNS


%% remove FRBs
FF= TNS.type~="FRB";
TNS = TNS(FF,:);
size(TNS,1)
% left with 195,747

%% save table
save -v7.3 TNS_NoFRB.mat TNS

%% add redshifts
R=imUtil.cat.match2Galaxies(TNS.ra,TNS.declination);

TNS = tools.table.addStructFields2TableCols(TNS, R);

Z = tools.array.selectFirstNotNaN(TNS.redshift, TNS.Z_PGC, TNS.Z_GLADE);

TNS.Z = Z;

%% save table
save -v7.3 TNS_withZ.mat TNS


  
%%

load TNS_withZ.mat;



% write table of TNS with redshifts
AllWithZ  = TNS(:,{'objid','name','type','ra','declination','redshift','Z'});
tools.table.sprintf_table(AllWithZ(1:5,:),'Format',{'%7d','%-7s','%-9s','%10.6f','%10.6f','%6.4f','%6.4f'},'IsLatex',true)
TableStr=tools.table.sprintf_table(AllWithZ,'Format',{'%7d','%-7s','%-9s','%10.6f','%10.6f','%6.4f','%6.4f'},'IsLatex',false);
FID = fopen('Table_TNS_WithZ.txt','w');
fprintf(FID,'%s',TableStr);
fclose(FID);

%%

% Skip this, since we are interested in all events

FlagZ=~isnan(TNS.Z); % <0.1; %06;    
TNSz=TNS(FlagZ,:);

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



%%
save -v7.3 Matched.mat MatchedRadio MatchedDist

sum(MatchedDist(:,1)<1)

    
%%
% select radio point sources
TNS.VLASS_Ftot = MatchedRadio(:,5);
TNS.VLASS_Fpeak= MatchedRadio(:,7);
TNS.VLASS_ErrFpeak = MatchedRadio(:,8);
TNS.VLASS_Dist     = MatchedDist(:,1);
Fps = ~isnan(MatchedRadio(:,5)) & MatchedRadio(:,5)< (MatchedRadio(:,7) + 3.* sqrt(MatchedRadio(:,6).^2 + MatchedRadio(:,8).^2));
% old compactness
Compactness = MatchedRadio(:,5) ./ (MatchedRadio(:,7) + 3.* sqrt(MatchedRadio(:,6).^2 + MatchedRadio(:,8).^2)); % compact are <1

% new compactness (C)

% for all VLASS
F_total = VLASS.Table.Ftot;
F_peak  = VLASS.Table.Fpeak;
DeltaF_total = VLASS.Table.ErrFtot;
DeltaF_peak  = VLASS.Table.ErrFpeak;

CV = ((F_total./F_peak) - 1) ./ ...
    (abs(F_total./F_peak) .* sqrt((DeltaF_total./F_total).^2 + (DeltaF_peak./F_peak).^2));



% for TNS
F_total = TNS.VLASS_Ftot;
F_peak  = TNS.VLASS_Fpeak;
DeltaF_total = MatchedRadio(:,6);
DeltaF_peak  = TNS.VLASS_ErrFpeak;


C = ((F_total./F_peak) - 1) ./ ...
    (abs(F_total./F_peak) .* sqrt((DeltaF_total./F_total).^2 + (DeltaF_peak./F_peak).^2));

TNS.Compactness = C;




% select clear associations (<1"):
Fd=MatchedDist(:,1)<1; 


%VLASS = cats.radio.VLASS1;
%FPS=~isnan(VLASS.Table.Ftot) & VLASS.Table.Ftot<(VLASS.Table.Fpeak + 3.*sqrt(VLASS.Table.ErrFtot.^2 + VLASS.Table.ErrFpeak));
%sum(FPS)./size(VLASS.Catalog,1)
sum(CV<2)./size(VLASS.Catalog,1)
% 41% are point sources 
sum(C<2 & Fd)./sum(Fd)
% 37%

% redshift
Fz = TNS.Z<0.1;
% 57,346 matches

% discovered prior to 2017 (VLASS epoch 1):
TNS.discoverydate_DateTime =   datetime(TNS.discoverydate);

Ftime = TNS.discoverydate_DateTime.Year<2017;
% 12,700 matches 

%Fprs = Fps & Fd & Fz & Ftime;
%sum(Fprs)
% 10 matches

FF = C<2 & Fd & Fz & Ftime;
sum(FF)
% 15 matches

Fprs = MatchedDist(:,1)<1.0 & Ftime;
sum(Fprs)
% 95

Fprs = MatchedDist(:,1)<3.3 & Ftime;
sum(Fprs)
% 176

Fprs = MatchedDist(:,1)<3.3 & Ftime & C<2;
sum(Fprs)
% 68



% all are Ic-BL !
%TNS.type(FF)          
%TNS.name(FF)
% old (all Ic-BL)
%    "2016gox"
%    "2016coi"

% abs mag:
TNS.disc_abs_mag = TNS.discoverymag - (5.*log10(astro.cosmo.lum_dist(TNS.Z))-5);

%%

save -v7.3 TNS_Selected.mat TNS

%%

%Selected = TNS(Fprs,{'name','type','ra','declination','Z','discoverymag','internal_names','disc_abs_mag'});
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
sum(TNS.discoverydate_DateTime.Year==2016)./sum(TNS.discoverydate_DateTime.Year<2017)

% z distribution of all TNS prior to 2017 and have z
F=TNS.discoverydate_DateTime.Year<2017;  
%TNS=TNS(F,:);

Zed=(0:0.01:0.5);
[Nz]=histcounts(TNS.Z(F),Zed);       
Xz=(Zed(1:end-1)+Zed(2:end))./2;
bar(Xz,Nz) 
H=xlabel('z');
H.FontSize = 18;
H.Interpreter = 'latex';
H=xlabel('N');
H.FontSize = 18;
H.Interpreter = 'latex';



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


%% plot distance distribution

EdgesDist = (0:0.25:10);
DistCen   = (EdgesDist(1:end-1) + EdgesDist(2:end)).*0.5;
Ndist     = histcounts(MatchedDist(:,1),EdgesDist);
bar(DistCen, Ndist)
H = xlabel('Transients-Radio seperation [arcsec]');
H.FontSize = 18;
H.Interpreter = 'latex';
H = ylabel('Number');
H.FontSize = 18;
H.Interpreter = 'latex';

print TransRadioSep.eps -depsc2

EdgesDist = (0:0.25:10);
DistCen   = (EdgesDist(1:end-1) + EdgesDist(2:end)).*0.5;
Ndist     = histcounts(MatchedDist(Ftime,1),EdgesDist);
bar(DistCen, Ndist)
H = xlabel('Transients-Radio seperation [arcsec]');
H.FontSize = 18;
H.Interpreter = 'latex';
H = ylabel('Number');
H.FontSize = 18;
H.Interpreter = 'latex';


print TransRadioSep_Prior2017.eps -depsc2

EdgesCompact = (-1:0.1:1);
Ncompact     = histcounts(log10(Compactness), EdgesCompact);
CompactCen   = (EdgesCompact(1:end-1) + EdgesCompact(2:end)).*0.5;
bar(CompactCen, Ncompact)
H = xlabel('log$_{10}$ radio compactness');
H.FontSize = 18;
H.Interpreter = 'latex';
H = ylabel('Number');
H.FontSize = 18;
H.Interpreter = 'latex';

print RadioCompactness_All.eps -depsc2

EdgesCompact = (-1:0.1:1);
Ncompact     = histcounts(log10(Compactness(FF)), EdgesCompact);
CompactCen   = (EdgesCompact(1:end-1) + EdgesCompact(2:end)).*0.5;
bar(CompactCen, Ncompact)
H = xlabel('log$_{10}$ radio compactness');
H.FontSize = 18;
H.Interpreter = 'latex';
H = ylabel('Number');
H.FontSize = 18;
H.Interpreter = 'latex';

print RadioCompactness_Tran.eps -depsc2

CC = C(Fd);
EdgesC = (-3:0.5:8);
BinC   = (EdgesC(1:end-1) + EdgesC(2:end)).*0.5;
Nc=histcounts(CC,EdgesC);
Ncv=histcounts(CV,EdgesC);
bar(BinC, [Ncv(:)./sum(Ncv), Nc(:)./sum(Nc)]);
H = xlabel('Compactness');
H.FontSize = 18;
H.Interpreter = 'latex';
H = ylabel('Fraction');
H.FontSize = 18;
H.Interpreter = 'latex';

legend('All','TNS');

print RadioCompactness_All.eps -depsc2


%%
TNS_FF   = TNS(C<2 & MatchedDist(:,1)<3.3 & Ftime, :);
TNS_FF   = sortrows(TNS_FF, {'VLASS_Dist'});

% match to external catalogs:
NVSS  = cats.radio.NVSS;
FIRST = cats.radio.FIRST;
EROSITA = cats.X.eRosita;
Nff     = size(TNS_FF,1);
TNS_FF.NVSS = nan(Nff,1);
TNS_FF.FIRST = nan(Nff,1);
TNS_FF.EROSITRA = nan(Nff,1);

for Iff=1:1:Nff
    DistAll = celestial.coo.sphere_dist_fast(TNS_FF.ra(Iff)./RAD, TNS_FF.declination(Iff)./RAD, NVSS.Cat(:,1), NVSS.Cat(:,2));
    [MinDist, MinInd] = min(DistAll);
    MinDist = MinDist.*RAD.*3600;
    if MinDist<10
        TNS_FF.NVSS(Iff) = NVSS.Cat(MinInd, 5);
    end
    
    DistAll = celestial.coo.sphere_dist_fast(TNS_FF.ra(Iff)./RAD, TNS_FF.declination(Iff)./RAD, FIRST.Cat(:,1), FIRST.Cat(:,2));
    [MinDist, MinInd] = min(DistAll);
    MinDist = MinDist.*RAD.*3600;
    if MinDist<1.5
        TNS_FF.FIRST(Iff) = FIRST.Cat(MinInd, 4);
    end
    
    DistAll = celestial.coo.sphere_dist_fast(TNS_FF.ra(Iff)./RAD, TNS_FF.declination(Iff)./RAD, EROSITA.Catalog(:,1), EROSITA.Catalog(:,2));
    [MinDist, MinInd] = min(DistAll);
    MinDist = MinDist.*RAD.*3600;
    if MinDist<7
        TNS_FF.EROSITA(Iff) = EROSITA.Catalog(MinInd, 14); % flux 0.2-2.3
    end

end

TNS_FF.EROSITA(TNS_FF.EROSITA == 0) = NaN;
TNS_FF.nuFnu = 3e9.*TNS_FF.VLASS_Fpeak.*4.*pi.*(astro.cosmo.lum_dist(TNS_FF.Z).*constant.pc).^2 .*1e-26;

% sort by no FIRST/NVSS counterparts:
F1 = TNS_FF.VLASS_Dist<=1;
F1_1 = isnan(TNS_FF.NVSS) & isnan(TNS_FF.FIRST) & F1;
F1_2 = (~isnan(TNS_FF.NVSS) | ~isnan(TNS_FF.FIRST)) & F1;
F2_1 = isnan(TNS_FF.NVSS) & isnan(TNS_FF.FIRST) & ~F1;
F2_2 = (~isnan(TNS_FF.NVSS) | ~isnan(TNS_FF.FIRST)) & ~F1;

TNS_FF1 = [TNS_FF(F1_1,:); TNS_FF(F1_2,:); TNS_FF(F2_1,:); TNS_FF(F2_2,:)];

%%

Selected = TNS_FF1(:,{'name','type','ra','declination','Z','discoverymag','disc_abs_mag', 'VLASS_Fpeak', 'VLASS_ErrFpeak', 'VLASS_Dist', 'nuFnu', 'NVSS', 'FIRST','EROSITA'});
tools.table.sprintf_table(Selected,'Format',{'%-7s','%-9s','%10.6f','%10.6f','%6.4f','%4.1f','%6.1f', '%6.1f', '%6.1f', '%3.1f', '%5.1e', '%5.1f', '%5.1f', '%4.1e'},'IsLatex',true, 'NoData',true)

