function x = FleshOut(x,lx,ASelect,u)
% Determine unrhymed syllables, given a framework of rhymed syllables

T = length(x);
fTQ = [0, find(x), T+1];
for k=1:(length(fTQ)-1)
    if fTQ(k+1) > fTQ(k) + 1
        Suffices = [];
        if k==length(fTQ)-1
            suffix = ones(lx.N,1);
        else
            suffix = zeros(lx.N,1);  suffix(x(fTQ(k+1))) = 1;
        end
        for t=(fTQ(k+1)-1):-1:(fTQ(k)+1)
            suffix = ASelect(u(t),:)' .* (lx.M * suffix);
            Suffices = [Suffices, suffix];  %  note, reverse order
        end
        for t=(fTQ(k)+1):(fTQ(k+1)-1)
            if t==1
                x(t) = SampFrom(cumsum(Suffices(:,end)'));
            else
                x(t) = SampFrom(cumsum(lx.M(x(t-1),:) .* Suffices(:, end)'));
            end
            Suffices(:,end)=[];
        end
    end
end