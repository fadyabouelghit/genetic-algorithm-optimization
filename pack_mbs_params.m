function [mbs_params, antennaObjectMbs, containsMbs, numBS] = pack_mbs_params(xs, ys, mbs_height, mbs_power, mbsAntenna)
% xs, ys: column vectors of equal length
numBS       = numel(xs);
containsMbs = double(numBS > 0);
mbs_x       = xs(:).';
mbs_y       = ys(:).';
mbs_height  = mbs_height * ones(1, numBS);
mbs_power   = mbs_power  * ones(1, numBS);
mbs_params  = [mbs_x; mbs_y; mbs_height; mbs_power];
antennaObjectMbs = repmat(mbsAntenna, 1, numBS);  % replicate your template
end
