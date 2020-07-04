function SpellOut(x,lx,rhysch)
% Input syllable codes and output a formatted verse

for t=1:length(x)
    if lx.isLastSyl(x(t))
        fprintf('%s ', lx.word{lx.sylToWord(x(t))});
        if any(rhysch.eols == t)
            fprintf('\n');
        end
    end
end
fprintf('\n');