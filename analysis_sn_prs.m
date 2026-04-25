%%

% downloaded on 6-Mar-2026
[TNS]=VO.TNS.downloadAll
size(TNS,1)
% 190,894

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

% Skip this, since we are interested in all events

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



%%
save -v7.3 Matched.mat MatchedRadio MatchedDist

sum(MatchedDist(:,1)<1)
% 1342 unique matches (with Z) ie 2420 (w and w/o Z)

    
%%
% select radio point sources
TNS.VLASS_Ftot = MatchedRadio(:,5);
TNS.VLASS_Fpeak= MatchedRadio(:,7);
TNS.VLASS_ErrFpeak = MatchedRadio(:,8);
TNS.VLASS_Dist     = MatchedDist(:,1);
Fps = ~isnan(MatchedRadio(:,5)) & MatchedRadio(:,5)< (MatchedRadio(:,7) + 3.* sqrt(MatchedRadio(:,6).^2 + MatchedRadio(:,8).^2));
% old compactness
Compactness = MatchedRadio(:,5) ./ (MatchedRadio(:,7) + 3.* sqrt(MatchedRadio(:,6).^2 + MatchedRadio(:,8).^2)); % compact are <1

% new compactness

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


% 2236 selected (with Z) / 2912 (all) 


% select clear associations (<1"):
Fd=MatchedDist(:,1)<1; 
% 1343 unique matches (with Z) / 2420 (all)


%VLASS = cats.radio.VLASS1;
%FPS=~isnan(VLASS.Table.Ftot) & VLASS.Table.Ftot<(VLASS.Table.Fpeak + 3.*sqrt(VLASS.Table.ErrFtot.^2 + VLASS.Table.ErrFpeak));
%sum(FPS)./size(VLASS.Catalog,1)
sum(CV<2)./size(VLASS.Catalog,1)
% 41% are point sources (with Z) / 81% (all)
sum(C<2 & Fd)./sum(Fd)
% 37%

% redshift
Fz = TNS.Z<0.1;
% 53,266 matches

% discovered prior to 2017 (VLASS epoch 1):
Ftime = TNS.discoverydate.Year<2017;
% 3973 matches (with Z) / 7588 (all)

FF = C<2 & Fd & Fz & Ftime;
sum(FF)
% 10 matches

% we are interested in all events
FF = Fd & Ftime;           
sum(FF)
% 86 matches (with 1.25'') / 120 with 3.3''
s

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

Selected = TNS_FF(:,{'name','type','ra','declination','Z','discoverymag','disc_abs_mag', 'VLASS_Fpeak', 'VLASS_ErrFpeak', 'VLASS_Dist', 'nuFnu', 'NVSS', 'FIRST','EROSITA'});
tools.table.sprintf_table(Selected,'Format',{'%-7s','%-9s','%10.6f','%10.6f','%6.4f','%4.1f','%6.1f', '%6.1f', '%6.1f', '%3.1f', '%5.1e', '%5.1f', '%5.1f', '%4.1e'},'IsLatex',true, 'NoData',true)

