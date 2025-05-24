function params = mergeParams(defaults, user)
    params = defaults;
    userFields = fieldnames(user);
    for i = 1:length(userFields)
        params.(userFields{i}) = user.(userFields{i});
    end
end
