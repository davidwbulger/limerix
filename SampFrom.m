function lix = SampFrom(cp)
% Use bisection search to sample from the given cumulative probability vector.
% DOES assume all probs are nonnegative. DOESN'T assume they are normalised.
L = length(cp);
p = rand*cp(L);
lix = 1;  hix = L;
while lix<hix
    mix = floor((lix+hix)/2);
    if cp(mix)<p
        lix = mix+1;
    else
        hix = mix;
    end
end