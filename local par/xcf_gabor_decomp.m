function mp = xcf_gabor_decomp(data,varargin)

%% define parameter
% data:   date for MP. Should be in the follow matrix: channelNum/trialLen/trialNum;
fieldlist = { 'srate'        'integer' []   [];           % sampling rate
              'mpIteration'  'integer' []   50;           % iteration number for MP
              'mpName'       'string' []    'MPtest';     % name of the dataset
              'mpPath'       'string'  []   cd;           % folder path to save the results
              };   
g = xcf_finputcheck( varargin, fieldlist);
if ischar(g), error(g); end
%xcf_newfolder(g.mpPath);
nbchan  = size(data,1);
nbpnt   = size(data,2);
nbtrial = size(data,3);
if mod(log2(nbpnt),1)~=0; error('len should be power of 2'); end

%%
% only perform the MP when the MP folder is not exist
if ~exist(platformSpecificNameMPP([g.mpPath '/' g.mpName '/GaborMP/mp0.bok.000']),'file')
    disp('Running MP decomposition...');
    dataMP   = nan(nbpnt,nbtrial*nbchan);
    for ch = 1:nbchan
        dataMP(:,nbtrial*(ch-1)+1:nbtrial*ch) = double(squeeze(data(ch,:,:)));
    end
    
    % Perform Gabor decomposition
    importData(dataMP,[g.mpPath '/'],[g.mpName '/'],[1 nbpnt],g.srate);       % import Data
    runGabor([g.mpPath '/'],[g.mpName '/'],nbpnt,g.mpIteration);            % Perform Gabor decomposition
    
    % retrieve information
    mp   = [];
    mp.g = g;
    gbtmp = getGaborData([g.mpPath '/'], g.mpName, 1); % Retrieve information
    for ch = 1:nbchan
        mp.gb(ch,:) = gbtmp(nbtrial*(ch-1)+1:nbtrial*ch);
    end
    
    % save the data
    mp.fs     = g.srate;
    mp.L      = nbpnt;
    save([g.mpPath '/' g.mpName '/mp.mat'],'mp');
    disp('calMP done!');
else
    disp('MP decomposition data already exists. Skipping MP decomposition...');
end