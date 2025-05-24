function clamped = clampToBounds(values, bounds)
    clamped = min(max(values, bounds(:,1)'), bounds(:,2)');
end
