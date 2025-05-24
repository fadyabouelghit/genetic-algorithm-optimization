function [x, y] = maternPP(areaSize, numPoints, minDist, maxAttempts)
    % Generate repulsive point pattern
    
    points = zeros(numPoints, 2);
    attempts = 0;
    currentPoint = 1;
    
    while currentPoint <= numPoints && attempts < maxAttempts
        newPoint = rand(1,2) .* areaSize;
        

        % use if biasing initial population is required. 

        % biasFactor = 1;  % Higher means stronger bias toward top-left
        % biasedX = rand^biasFactor * areaSize(1);
        % biasedY = (1 - rand^biasFactor) * areaSize(2);
        % newPoint = [biasedX, biasedY];

        if currentPoint == 1
            points(currentPoint,:) = newPoint;
            currentPoint = currentPoint + 1;
        else
            distances = pdist2(newPoint, points(1:currentPoint-1,:));
            if all(distances > minDist)
                points(currentPoint,:) = newPoint;
                currentPoint = currentPoint + 1;
            end
        end
        attempts = attempts + 1;
    end
    
    % Fill remaining points randomly if needed
    while currentPoint <= numPoints
        points(currentPoint,:) = rand(1,2) .* areaSize;
        currentPoint = currentPoint + 1;
    end
    
    x = points(:,1);
    y = points(:,2);
end
