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

FF = Fps & Fd & Fz & Ftime;
sum(FF)
% 10 matches

% all are Ic-BL !
%TNS.type(FF)          
%TNS.name(FF)
% old (all Ic-BL)
%    "2016gox"
%    "2016coi"

% abs mag:
TNS.disc_abs_mag = TNS.discoverymag - (5.*log10(astro.cosmo.lum_dist(TNS.Z))-5);

save -v7.3 TNS_Selected.mat TNS

%%
Selected = TNS(FF,{'name','type','ra','declination','Z','discoverymag','internal_names','disc_abs_mag'})
tools.table.sprintf_table(Selected,'Format',{'%-7s','%-9s','%10.6f','%10.6f','%6.4f','%4.1f','%-35s','%6.1f'},'IsLatex',true)

