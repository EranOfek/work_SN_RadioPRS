%%

% downloaded on 3-Mar-2026
[Tsn]=VO.TNS.downloadAll;

save -v7.3 Tsn.mat Tsn

% 190,740
  
%%

load Tsn.mat;

%%

FlagZ=Tsn.redshift<0.1; %06;    
Tsn=Tsn(FlagZ,:);

% 17,397

%%
RAD = 180./pi;

SearchRadius = 10;  % [arcsec]

Nsn = size(Tsn,1);
tic;
for Isn=1:1:Nsn
	  [Isn, Nsn]
    [CatRadio, ColRadio, ~, Dist] = catsHTM.cone_search('VLASSep1', Tsn.ra(Isn)./RAD, Tsn.declination(Isn)./RAD, SearchRadius);
     
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

    % 773 matches, 677 unique matches 

    
%%
% select radio point sources
Fps = ~isnan(MatchedRadio(:,5)) & MatchedRadio(:,5)./MatchedRadio(:,7)<(1+3.* MatchedRadio(:,6));

% 511 left

% select clear associations (<1"):
Fd=MatchedDist(:,1)<1; 

% redshift
Fz = Tsn.redshift<0.1;

% discovered prior to 2017 (VLASS epoch 1):
Ftime = Tsn.discoverydate.Year<2017;

FF = Fps & Fd & Fz & Ftime;
sum(FF)

% all are Ic-BL !
Tsn.type(FF)          
Tsn.name(FF)
%    "2016gox"
%    "2016coi"

