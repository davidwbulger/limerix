function rhysch = MakeRhymeScheme(infile)
% Reads in a human-readable rhyming scheme description, and outputs a
% structure representing it in the form required by OutputRhyme.
% For instance, if infile contains
% .-..-..A
% .-..-..A
% .-..B
% .-..B
% .-..-..A
% then this function should output struct('stress', [3 2 1 1 2 1 1 6 3 2 1 1 2 1 1 6 3 2 1 1 6 3 2 1 1 6 3 2 1 1 2 1 1 6], 'sylsets', {{[8, 16, 34], [21, 26]}}, 'eols', [8,16,21,26,34]);
% 'Stress' codes: 1=unstressed; 2=stressed; 3=first syl unst; 4=first syl str; 5=last syl unst; 6=last syl str; 7=monosyl unst; 8 = monosyl str.
% Use lowercase letters for rhymed unstressed syllables.

fid = fopen(infile, 'r');
firstsyl = false(1,0); lastsyl = firstsyl; glyphs = '';
l = fgetl(fid);
while ~isequal(l, -1),
    glyphs = [glyphs, l];
    firstsyl = [firstsyl, true, false(1, length(l)-1)];
    lastsyl = [lastsyl, false(1, length(l)-1), true];
    l = fgetl(fid);
end
isupper = (lower(glyphs)~=glyphs);
stressed = isupper | (glyphs=='-') | (glyphs=='_');
rimes = upper(glyphs);
sylsets = {};
while any(isletter(rimes)),
    rx=find(isletter(rimes),1);
    rX=find(rimes==rimes(rx));
    sylsets = {sylsets{:}, rX};
    rimes(rX) = '&';
end
rhysch = struct('stress', 1 + stressed + 2*firstsyl + 4*lastsyl, 'sylsets', {sylsets}, 'eols', find(lastsyl));