function value = reflectToBounds(val, bound)
% Reflect exceeded values back into bounds
    range = bound(:,2) - bound(:,1);

    over = val > bound(:,2)';
    under = val < bound(:,1)';

    val(over) = 2*bound(over,2)' - val(over);
    val(under) = 2*bound(under,1)' - val(under);

    needsCorrection = over | under;
    if any(needsCorrection)
        val(needsCorrection) = reflectToBounds(val(needsCorrection), bound(needsCorrection,:));
    end
    
    value = val;
end
