# limerix
Create random limericks based on source texts.

Look at the article (the .pdf file) to read an overview of the project and see some output samples.

To run this yourself in Matlab, try:

    load NorthangerAbbeyLexicon.mat
    rhsc = MakeRhymeScheme('limerick.txt');
    AvoidSelfRhyme(lexicon, rhsc, 4, 30);
  
In order to get MakeFrags.m working (i.e., to train the Markov chain with anything other than Northanger Abbey),
you'll need to get hold of the CMU pronunciation dictionary. There may be some other steps required;
get in touch (david.bulger@mq.edu.au) if you have trouble.

I'm just quoting the article's abstract:

Homogeneous Markov chains with discrete state-space and time are very straightforward to
simulate, due to the memorylessness property. In this work, we consider sampling from the normalised
restriction of a Markov chains joint distribution to a subset defined by two kinds of constraints:
univariate constraints, restricting individual chain variables to belong to certain sets, and multivariate
constraints, restricting non-adjacent groups of chain variables to be unequal but equivalent under a
predefined equivalence relation. The multivariate constraints, in particular, present some challenges. An
efficient sampling method involving Metropolis-Hastings with some pre-tabulated distributions is described.
We demonstrate an application to text generation, and in particular the random generation of rhyming, scanning
lyrics, as a component of an algorithmic songwriting project. In this setting, the states of the Markov chain
are syllables (that know which words they are from, so that for instance the final syllables of perilous and
marvellous are two different states). The univariate constraints enforce the metre of the text, by indicating
which syllables must be stressed or unstressed (for correct scansion) and which syllables must be word-final
(so that words do not straddle lines). The multivariate constraints enforce a rhyming scheme, by requiring
certain groups of syllables to rhyme without being equal.

This work relies on the Carnegie Mellon Pronouncing Dictionary (CMUdict; Weide, 1998) and sample text.
CMUdict contains over 130000 entries, comprising words, proper nouns et cetera as used in spoken American
English. This dictionary indicates the phonemic and stress pronunciation of each entry.
A Markov chain including the entire dictionary would require over 330000 syllable states. We have used
sample texts, both to reduce the state space, and to build a first-order Markov m odel of word sequence, or
equivalently, syllable sequence. The problem of interest is to sample from the Markov chain, subject to the
constraints imposed by a given rhyming scheme.

Our algorithm firstly samples the multivariately constrained chain variables. After location of a feasible
subsequence of rhymed syllables, ratios of joint probabilities are calculated and compared in a Metropolis-Hastings
algorithm, where each candidate arises from the incumbent by replacing one set of rhyming syllables with
another drawn from a pre-tabulated distribution. After a burn-in period, the skeleton provided by the rhymed
syllables is fleshed out by sampling the unrhymed syllables according to the Markov chain, subject to metrical
and boundary conditions.
