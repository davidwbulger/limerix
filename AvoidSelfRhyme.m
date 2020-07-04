function xset = AvoidSelfRhyme(lx, rhysch, numverses, cheatafter)
% Output random text matching the rhyming scheme 'rhysch', based on stats in lx.
% 'Stress' codes: 1=unstressed; 2=stressed; 3=first syl unst; 4=first syl str; 5=last syl unst; 6=last syl str; 7=monosyl unst; 8 = monosyl str. (Note: these are NOT exclusive.)
% Also, numverses is the number of outputs required (with a new output
% whenever each rime has changed) and cheatafter is how many proposals in a row
% it will reject before accepting automatically.
% Therefore, call thus:
%   Do either
%     MakeFrags;  %  creates the lexicon (24s) (just over 10min for MobyDick.txt)
%   or
%     load MobyDickLexicon.mat;  %  (0.4s)
%   rhsc = MakeRhymeScheme('limerick.txt');  %  (0.03s)
%   AvoidSelfRhyme(lexicon, rhsc, 10, 30);  %  (~10s)

% Note:
%   This is based on the earlier 'OutputRhyme', but rebuilt so as to avoid
%   self-rhymes.
%
%   We expect the problem of improper scansion of 1-syl words to persist.
%   (Formally, the single syllable in a one-syllable word is 'stressed',
%   but in natural speech, some singletons must be or not be stressed
%   depending on grammar and semantics; this process ignores that.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function pr = PropRat
        % Proposal density ratio. The proposal density is not symmetric, so
        % the Metropolis-Hastings algorithm requires the ratio of the
        % proposal densities between incumbent and candidate.
        % Because the proposal density is just the normalised restriction
        % of mu^K to the feasible sequences, the proposal ratio (assuming
        % both sequences are feasible) is simply the ratio of mu^K.
        ks = length(cand);
        incumbent = x(Q(r,:));
        pr = lx.mu(incumbent(1))/lx.mu(cand(1));
        for kx = 2:ks
            pr = pr * lx.mu(incumbent(kx))/lx.mu(cand(kx));
        end
    end

    function lr = LikeRat(rho, sigma)
        % Likelihood ratio between candidate (in longcand) and incumbent (in x).
        if x(rho)
            incprovec = zeros(1,lx.N);
            canprovec = zeros(1,lx.N);
            incprovec(x(rho)) = lx.mu(x(rho));
            canprovec(longcand(rho)) = lx.mu(longcand(rho));
        else
            incprovec = lx.mu .* ASelect(u(rho),:);
            canprovec = incprovec;
        end
        for tx=(rho+1):sigma
            incprovec = incprovec * lx.M;
            canprovec = canprovec * lx.M;
            if tx<sigma
                if TQ(tx)
                    incprovec(1:(x(tx)-1)) = 0;
                    incprovec((x(tx)+1):end) = 0;
                    canprovec(1:(longcand(tx)-1)) = 0;
                    canprovec((longcand(tx)+1):end) = 0;
                else
                    incprovec = incprovec .* ASelect(u(tx),:);
                    canprovec = canprovec .* ASelect(u(tx),:);
                end
            end
        end
        lr = canprovec(longcand(sigma)) / incprovec(x(sigma));
    end

    function s = JointSamp(Kx, ux)
        % This returns a vector of K syllables distributed according to the
        % normalised restriction of mu^K to those vectors that rhyme and are
        % distinct. It relies on h, sylvie and vprobs as set up in the main
        % function, so possibly it needs to be 'nested'.
        
        % Firstly choose the rime:
        vx = SampFrom(vprobs{ux}(:,Kx));
        
        % Now iteratively choose the highest, second highest &c.:
        s = [zeros(Kx-1,1); SampFrom(h{vx,ux}(:, Kx))];
        for kx = (Kx-1):-1:1
            s(kx) = SampFrom(h{vx,ux}(1:(s(kx+1)-1), kx));
        end
        
        % Now randomly permute:
        s = sylvie{vx,ux}(s(randperm(length(s))));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING AND INITIALISATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rhyming scheme is initially specified in this form:
% rhysch.stress is an integer vector indicating which syllables should be stressed, which must be word-initial, and which must be word-terminal;
% rhysch.sylsets is a cell array of vectors of syllable indices ('sylsets') which must rhyme;
% rhysch.eols determines the location of line breaks in the final output.
verbose = false;

u = rhysch.stress(:);  %  which set A_u_t each X_t has to belong to, according to univariate constraints.
T = length(u);
N = lx.N;  %  number of states, i.e., syllables in dictionary
U = 8;  %  number of 'stress' codes
V = max(lx.rime);  %  number of rimes in lexicon
R = length(rhysch.sylsets);
Q = zeros(R,T);
K = zeros(R,1);  %  number of syllables involved in each rhyme constraint
for r=1:R
    Q(r,rhysch.sylsets{r}) = 1;
    K(r) = length(rhysch.sylsets{r});
end
Q = logical(Q);
maxK = max(K);
TQ = sum(Q);
if any(TQ>1)
    error(['Multiple multivariate constraints at X_', num2str(find(TQ>1,1))]);
else
    TQ = logical(TQ);
end

% Initialise utility arrays for matrix calculations

% Note that here we're allowing any one-syllable word to appear in an
% unstressed position, as well as a stressed position. In fact, it
% shouldn't be that simple. There should be three types of 1-syl words,
% allowed in stressed or unstressed positions or either, and we don't have
% that info. Maybe write to CMU people?

% Anyway, this is a binary matrix with eight rows, where the jth row just
% has a '1' for each syllable in the jth stress category:
ASelect = [~lx.stress | (lx.isFirstSyl & lx.isLastSyl); lx.stress; ~lx.stress & lx.isFirstSyl; lx.stress & lx.isFirstSyl; ~lx.stress & lx.isLastSyl; lx.stress & lx.isLastSyl; ~lx.stress & lx.isFirstSyl & lx.isLastSyl; lx.stress & lx.isFirstSyl & lx.isLastSyl];
%ASelect = [~lx.stress ; lx.stress; ~lx.stress & lx.isFirstSyl; lx.stress & lx.isFirstSyl; ~lx.stress & lx.isLastSyl; lx.stress & lx.isLastSyl; ~lx.stress & lx.isFirstSyl & lx.isLastSyl; lx.stress & lx.isFirstSyl & lx.isLastSyl];

% For each rime v, the vth row is a binary vector of all the rhyming syllables:
BSelect = zeros(V,N);
for n=1:N
    BSelect(lx.rime(n), n) = 1;
end

% Determine the fixed boundary points for each set of rhymes:
RhoSigma = zeros(2,R);
for r = 1:R
    TQr = sort(rhysch.sylsets{r});  %  probably already sorted but be safe
    TQs = find(TQ-Q(r,:));
    rho = [1, TQs];  rho(rho>min(TQr)) = [];  rho = max(rho);
    sigma = [TQs, T];  sigma(sigma<max(TQr)) = [];  sigma = min(sigma);
    RhoSigma(:,r) = [rho; sigma];
end
% The rth column is just the endpoints of the span containing the
% syllables of the rth ryme group. The lower value is the last fixed
% syllable BEFORE the rhyme group starts, or 1 if this rhyme is the
% earliest in the scheme. Similarly, the higher value (RhoSigma(2,r)) is
% the first fixed value AFTER it ENDS, or T if this is the last rhyme
% appearing in the scheme.

% Determine rough-estimate distributions on rimes. For each positive
% integer k up to maxK, and for each stress category ux, we want a
% probability distribution on the rimes, as they arise iid ~mu, conditional
% on the states all differing but rhyming and having stress category ux. In
% practice, most stress categories will not be needed (e.g., for limericks,
% only category 6 is used) and so we do not compute them all. Note that we
% are NOT allowing for cases where one rhyming group in the scheme includes
% multiple stress groups (as might occur if internal rhymes are schemed).
%h = cell(U, maxK);
% DO THE STRESS CODES AS A CELL ARRAY RATHER THAN A 3RD DIMENSION, SINCE WE
% PROBABLY DON'T NEED MOST OF THEM.
% ALSO, STORING OUT TO MAX # RHYMES PER RIME FOR EACH RIME WOULD BE WASTEFUL
h = cell(V, U);  %  initially, this will be the unnormalised PMF. V=#rimes in lex, maxK=max#matched syls in scheme, U=#stress codes
% Note that we don't bother with the factorials, since the probabilities
% are unnormalised, and for fixed k the effect will be a constant ratio.
sylvie = cell(V,U);  %  syllables rhyming with v (and with stress code u)
vprobs = cell(U,1);  %  for each stress code and each k, this gives the (unnormalised) CDF of the rime, v
for ux = unique(u(TQ))  %  that is, for each stress group appearing as a rhymed syllable  --  was unique(u(find(TQ)))
    vprobs{ux} = zeros(V, maxK);
    for v = 1:V
        sylvie{v,ux} = find((lx.rime'==v) & ASelect(ux,:));
        if length(sylvie{v,ux})
            p = lx.mu(sylvie{v,ux})';  %  stationary probabilities of matching syllables, dictionary order
            h{v,ux} = [cumsum(p), zeros(length(p), maxK-1)];
            for k = 2:maxK
                h{v,ux}(k:end,k) = cumsum(p(k:end) .* h{v,ux}((k-1):(end-1), k-1));
            end
            vprobs{ux}(v,:) = h{v,ux}(end,:);
        else
            vprobs{ux}(v,:) = 0;
        end
    end
    vprobs{ux} = cumsum(vprobs{ux});
end
% Now the columns of h are unnormalised CDFs, but proportional for each k.

% Initialise multivariately constrained chain vars:
x = zeros(1,T);
for r=1:R
    prob = 0;
    while ~prob
        x(Q(r,:)) = JointSamp(K(r), unique(u(Q(r,:))));
        % The main idea here is to ensure that the word that the
        % rhyming syllable belongs to actually fits into the rhyming
        % scheme. But the lexicon doesn't easily give us the start &
        % end position of the current word, and anyway, there must
        % certainly be words in the lexicon that appear only once, or
        % at least only ever preceded/followed by a certain other word.
        % So instead the plan is, after fitting each rhyme, we'll
        % verify that the overall probability is still positive, and if not,
        % we'll fit that rhyme again. It is conceivable that one rhyme
        % which has positive probability under the metric scheme would
        % nevertheless not allow any choice for the subsequent rhyme,
        % but that seems unlikely unless the rhyming scheme is very
        % cramped.
        fprintf('Trying to initialise rhyme #%d...\n', r);
        if x(1)
            provec = zeros(1,lx.N);
            provec(x(1)) = 1;
        else
            provec = ASelect(u(RhoSigma(1,r)),:);  %  Note that magnitudes are irrelevant here
        end
        for t=2:find(x,1,'last')
            provec = (provec * lx.M) .* ASelect(u(t),:);
            if x(t)
                provec(1:(x(t)-1)) = 0;
                provec((x(t)+1):end) = 0;
                prob = provec(x(t));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE MAIN LOOP TO SAMPLE THE RHYMED SYLLABLES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Should just need to change this section ... although we may want to add a
% section above to precalculate the n-step matrices for the stretches
% between rhymed syllables. (Or is that already done above?)
if verbose, fprintf('Entering Metropolis sampling phase.\n'); end
% We run Metropolis until each rime has changed at least once. The
% following vector keeps track of this stopping criterion:
versecount = 0;
xset = zeros(numverses, T);
while versecount < numverses
    changed = false(R,1);
    cheatcount = 0;
    while any(~changed)
        r = ceil(R*rand);
        if verbose, fprintf('Proposed change to rhyme #%d (out of %d): ', r, R); end
        
        ux = unique(u(Q(r,:)));
        if length(ux)>1
            error('Cross stress rhymes in scheme; rethink everything!');
        end
        cand = JointSamp(K(r), ux);
        longcand = x;
        longcand(Q(r,:)) = cand;
        if cheatcount >= cheatafter || rand < LikeRat(RhoSigma(1,r),RhoSigma(2,r))*PropRat  %  LikeRat(x,cand,r)*PropRat(x,cand,r)
            if cheatcount >= cheatafter
                fprintf('cheating\n');
            end
            ThisChanged = true;
            x = longcand;
            cheatcount = 0;
        else
            ThisChanged = false;
            cheatcount = cheatcount + 1;
        end
        
        changed(r) = changed(r) || ThisChanged;
        if verbose
            if ThisChanged
                fprintf('accepted. ');
            else
                fprintf('rejected. ');
            end
            fprintf('So far, %d rhymes changed.\n', sum(changed));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAMPLE THE UNRHYMED SYLLABLES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose, fprintf('Sampling the unrhymed syllables.\n\n'); end
    xo = FleshOut(x, lx, ASelect, u);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FINALLY, OUTPUT THE GENERATED TEXT:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SpellOut(xo,lx,rhysch);
    versecount = versecount + 1;
    xset(versecount,:) = xo;
end
end