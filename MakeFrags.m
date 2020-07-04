% LOAD DICTIONARY
tic;

debugging = false;
if debugging
    dictfn = 'SubDict.txt';  %  Actually drastically abridged (for debugging)
    samplefn = 'ShortSample.txt';
    outfn = 'lexiconShort.mat';
else
    dictfn = 'DictReduced.txt';  %  Only reduced by removing header material
    opna =   'NorthangerAbbey';  %  'MobyDick';  %  'Jeeves';
    samplefn = ['TextSamples/', opna, '.txt'];
    outfn = [opna, 'Lexicon.mat'];
end

fid = fopen(dictfn, 'r');
dictionary = struct([]);
reli = fgetl(fid);  %  'REad LIne' = first line of dictionary
while ischar(reli)
    ennu = length(dictionary) + 1;  %  entry number
    toks = strsplit(reli);
    dictionary(ennu).word = toks{1};
    dictionary(ennu).rime={};
    dictionary(ennu).stress=[];
    for k=2:length(toks)
        if any('AEIOU'==toks{k}(1))
            dictionary(ennu).rime=[dictionary(ennu).rime,{toks{k}(1:2)}];
            dictionary(ennu).stress=[dictionary(ennu).stress, toks{k}(end)>'0'];  %  conflating primary ('1') & secondary ('2') stress
        elseif ~isempty(dictionary(ennu).rime)
            dictionary(ennu).rime{end}=[dictionary(ennu).rime{end},toks{k}];
        end
    end
    if ~mod(ennu,1000)
        fprintf('%d words read from the dictionary.\n', ennu);
    end
    reli = fgetl(fid);  %  next line of dictionary
end
fclose(fid);
toc,
fprintf('Loaded dictionary.\n');

% LOAD TEXT SAMPLE
% This assumes substantial preprocessing, so that each word is in capitals
% & on its own line. That's done in gvim; see Cleaning.txt.
fid = fopen(samplefn, 'r');
sample = cell(0,0);
reli = fgetl(fid);
while ischar(reli)
    sample = [sample, {reli}];
    reli = fgetl(fid);
    l = length(sample);
    if ~mod(l,1000)
        fprintf('%d words read from the sample.\n', l);
    end
end
fclose(fid);
toc,
fprintf('Loaded text sample.\n');

% DETERMINE USED WORDS & FRAGMENTS:
[Lia,Loc] = ismember(sample, {dictionary.word});
Loc = [Loc,0];  %  prob unnecessary, but standardises shape
[C,~,ic]=unique(Loc);
Loc=ic'-1;
dictionary = dictionary(C(2:end));  %  e.g., knowing there's at least one 0 allows this (discard words appearing in dictionary but not sample)
sylof = zeros(length(dictionary),1);  %  syllable offsets
for k=2:length(dictionary)
    sylof(k) = sylof(k-1) + length(dictionary(k-1).rime);  %  i.e., plus # of syllables in word
end
N=0; for k=1:length(dictionary), N=N+length(dictionary(k).rime); end  %  State count.
strime = cell(1,N);  %  state rime
ststr = zeros(1,N);  %  state stress
%stord = zeros(1,N);  %  state ordinal syllable in word
ifs = false(1,N);  %  is first syllable
ils = false(1,N);  %  is last syllable
stw = zeros(1,N);  %  map syllable (state ID) to word ID
for k=1:length(dictionary)
    p = 1 + sylof(k);
    l = length(dictionary(k).rime);
    %stord(p:(p+l-1)) = 1:l;
    ifs(p) = true;
    ils(p+l-1) = true;
    stw(p:(p+l-1)) = k;
    for j=1:length(dictionary(k).rime)
        strime(p:(p+l-1)) = dictionary(k).rime;
        ststr(p:(p+l-1)) = dictionary(k).stress;
    end
end
[~,~,rid] = unique(strime);  %  serial number for each unique rime
toc,
fprintf('Fragmented.\n');

% CREATE TRANSITION MATRIX:
stafreqy = 0; % observed state frequencies.
M = zeros(N);
for k=1:length(dictionary)
    l = length(dictionary(k).rime);
    if l > 1
        s = sylof(k);
        M((s+1):(s+l), (s+1):(s+l)) = diag(ones(l-1,1),1);  %  intralexical trans probs
    end
end
for k=2:length(Loc)
    if all(Loc([k-1,k]))  %  i.e., if both words are in the dictionary
        M(sylof(Loc(k-1))+length(dictionary(Loc(k-1)).rime),sylof(Loc(k))+1) = M(sylof(Loc(k-1))+length(dictionary(Loc(k-1)).rime),sylof(Loc(k))+1) + 1;
    end
end
% Now normalise the rows:
M = M ./ (sum(M,2) * ones(1,N));

% This leaves NaNs for words that are only seen at the ends of fragments.
% Their rows, I guess, can be constructed from observed word proportions
% (weighting only, of course, first syllables).
woco = hist(Loc(Loc>0), 1:max(Loc));  %  Word counts.
%for k=1:length(woco),  %  Table of words & counts
%    fprintf('%s  %d\n', dictionary(k).word, woco(k));
%end
woco = woco / sum(woco);
fisy = zeros(1,N);  %  observed proportions, first syllables
mu = zeros(1,N);  %  stationary distribution of M;
for k=1:length(woco)
    fisy(1+sylof(k)) = woco(k);
    mu((1:length(dictionary(k).rime))+sylof(k)) = woco(k);
end
oldmu = mu / sum(mu);
for j=1:N
    if isnan(M(j,1))
        M(j,:) = fisy;
    end
end
mu = oldmu * M;
mu = mu/sum(mu);
while mean(abs(mu-oldmu)) > 1e-6
    oldmu = mu;
    mu = oldmu * M;
    mu = mu/sum(mu);
end

toc,
fprintf('Transition matrix and stationary distribution determined.\n');

% Package it all up and save it:
lexicon = struct('word', {{dictionary.word}}, ... % the words, in standard orthography, used for output
    'N', N, ... % the number of states (one for each syllable of each word)
    'M', M, ... % the transition probability matrix
    'insta', fisy, ... % the initial state distribution; weighted on first syllables, based on word proportions in sample
    'mu', mu, ... % the stationary distribution, estimated
    'rime', rid, ... % the Rime ID of each syllable (by state number)
    'stress', ststr, ...
    'isFirstSyl', ifs, ...
    'isLastSyl', ils, ...
    'sylToWord', stw);
save(outfn, 'lexicon', '-v7.3');