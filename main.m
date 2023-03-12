% Tao Xie, Mar/11/2023
% contact: xie@neurotechcenter.org

addpath(genpath('MP_master'))   % Matching Pursuit method. For detail please see https://github.com/supratimray/MP
addpath('local par')

%% (1) Extracting slow-wave and higher frequency band power
% We directly load the pre-processed data with 25Hz sampling rate, please see the manuscript for detailed pre-processing procedure
% EEGh: broadband gamma envelope
% EEGl: slow-wave 
load('local par/pre-processed data.mat','EEGh','EEGl')

%% (2) Calculating the tau-modulation curve (TMC)
if false
    step  = 5;       % stepping every 5 sample points, which is equal to 0.2 second
    lag   = -31:32;  % width equal to 2.56 second
    L     = length(lag);
    loc   = L+1:step:EEGl.pnts-L-1;
    tmc      = [];
    tmc.tLag = lag./EEGh.srate;
    tmc.t    = EEGl.times(loc);
    rAll     = zeros(EEGh.nbchan,L,length(loc));
    hamDat   = hamming(L);
    for ch = 1:EEGl.nbchan
        ld = double(EEGl.data(ch,:));
        hd = double(EEGh.data(ch,:));
        r  = zeros(L,length(loc));
        for i = 1:length(loc)
            lSeg = ld(loc(i)+lag)'.*hamDat;
            for k = 1:L
                hSeg  = hd(loc(i)+lag+lag(k))'.*hamDat;
                r(k,i) = corr(lSeg,hSeg,'type','pearson');
            end
        end
        rAll(ch,:,:) = r;
    end
    tmc.r = rAll;
    save('local par/tmc.mat','tmc')
end

%% (3) Extracting the three features (strength, frequency, polarity) from the tau-modulation curve
load('local par/tmc.mat','tmc')
ftrTim = tmc.t;
sw     = 24; % step width
loc    = sw/2+1:length(ftrTim)-sw/2;

% extracting the modulation strength and polarity
ftrSNR = nan(EEGl.nbchan,length(tmc.t));
ftrPhs = nan(EEGl.nbchan,length(tmc.t));
rStr = [];
for i = 1:length(loc)
    for ch = 1:EEGl.nbchan
        dat = squeeze(tmc.r(ch,:,loc(i)-sw/2:loc(i)+sw/2));
        d   = detrend(squeeze(mean(dat,2)));
        ftrSNR(ch,loc(i)) = var(dat(:))/mean(var(dat,[],2)); % modulation strength
        ftrPhs(ch,loc(i)) = sign(mean(d(30:35)));            % modulation polarity
    end
end

% extracting the modulation frequency
% Matching Pursuit
if false
    dat  = zeros(size(tmc.r));
    for i = 1:length(loc)
        dat(:,:,loc(i)) = squeeze(mean(tmc.r(:,:,loc(i)-sw/2:loc(i)+sw/2),3));
    end
    ftrMPt = tmc.t(loc(1):25:loc(end));
    dat    = dat(:,:,loc(1):25:loc(end));
    xcf_gabor_decomp(dat,'srate',EEGl.srate,'mpIteration',5,'mpName','testData','mpPath',[cd '/MP']);
    
    % save the mp
    load([cd '/MP/testData/mp.mat'],'mp');
    mp.ftrMPt = ftrMPt;
    mp.fs     = EEGl.srate;
    mp.L      = size(dat,2);
    save('local par/mp.mat','mp');
end

% retrieve the ftrFre
load('local par/mp.mat','mp');
ftrFre = nan(EEGl.nbchan,length(mp.ftrMPt));
for ch = 1:EEGl.nbchan
    for i = 1:length(mp.ftrMPt)
        gbDat = mp.gb{ch,i}.gaborData;
        atom  = find(gbDat(2,:)~=0,1);
        if ~isempty(atom)
            ftrFre(ch,i)  = gbDat(2,atom).*(mp.fs/mp.L); % modulation frequency
        end
    end
end
ftrFret = mp.ftrMPt;

clearvars -except EEGh EEGl ftrFre ftrFret ftrTim ftrSNR ftrPhs mp

%% show the example
ch = 2; fsize = 16; xl = [10 50];
AnesInj   = EEGl.times(round(EEGl.event(3).latency))*1e-3/60; % AnestheticInjection
AnesStart = EEGl.times(round(EEGl.event(4).latency))*1e-3/60; % Anesthetized-Start
AnesEnd   = EEGl.times(round(EEGl.event(5).latency))*1e-3/60; % Anesthetized-End

figure('position',[10 10 700 800]);
% modulation strength
subplot(311);hold on;yl=[0.9 3];
plot(ftrTim*1e-3/60,ftrSNR(ch,:),'color',[.5 .5 .5],'linewidth',2); 
plot([1 1].*AnesInj,yl,'--k','linewidth',2);
plot([1 1].*AnesStart,[yl(1) yl(1)+diff(yl)/2],'--b','linewidth',1.5);
plot([1 1].*AnesEnd,  [yl(1) yl(1)+diff(yl)/2],'--b','linewidth',1.5);
xlabel('Time (min)'); ylabel('SNR'); title(['Monkey ' EEGl.subject ': ' EEGl.chanlocs(ch).labels])
xlim(xl); ylim(yl); set(gca,'fontsize',fsize,'xtick', [AnesInj AnesInj+10 AnesInj+20 AnesInj+30],'xticklabel',[0 10 20 30]); 

% modulation frequency
subplot(312);hold on;yl=[0 4];
plot(ftrFret*1e-3/60,ftrFre(ch,:),'color',[.5 .5 .5],'linewidth',2); 
plot([1 1].*AnesInj,yl,'--k','linewidth',2);
plot([1 1].*AnesStart,[yl(1) yl(1)+diff(yl)/2],'--b','linewidth',1.5);
plot([1 1].*AnesEnd,  [yl(1) yl(1)+diff(yl)/2],'--b','linewidth',1.5);
xlabel('Time (min)'); ylabel('Frequency (Hz)')
xlim(xl); ylim(yl); set(gca,'fontsize',fsize,'xtick', [AnesInj AnesInj+10 AnesInj+20 AnesInj+30],'xticklabel',[0 10 20 30]); 

% modulation polarity
subplot(313);hold on;yl=[-2 2];
phs = ftrPhs(ch,:);
plot(ftrTim*1e-3/60,phs,'color',[.5 .5 .5],'linewidth',1); 
plot(ftrTim(phs==1)*1e-3/60,phs(phs==1),'o','markerfacecolor','r','markeredgecolor','none'); 
plot(ftrTim(phs==-1)*1e-3/60,phs(phs==-1),'o','markerfacecolor','b','markeredgecolor','none'); 
plot([1 1].*AnesInj,yl,'--k','linewidth',2);
plot([1 1].*AnesStart,[yl(1) yl(1)+diff(yl)/2],'--b','linewidth',1.5);
plot([1 1].*AnesEnd,  [yl(1) yl(1)+diff(yl)/2],'--b','linewidth',1.5);
xlabel('Time (min)'); ylabel('Polarity')
xlim(xl); ylim(yl); set(gca,'fontsize',fsize,'xtick', [AnesInj AnesInj+10 AnesInj+20 AnesInj+30],'xticklabel',[0 10 20 30]); 

saveas(gca,'tau-modulation freatures','png')



