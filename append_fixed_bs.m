function [mbs_params, antennaObjectMbs, containsMbs, numMbs] = append_fixed_bs(mbs_params, antennaObjectMbs, extras, defaultAntenna)
% Append extra fixed (non-GA-controlled) base stations to an existing
% mbs_params / antennaObjectMbs pair produced by pack_mbs_params followed by
% the row-1<->row-2 swap used in optimize_base_station_ga.m.
%
% Inputs:
%   mbs_params       : 4×N current matrix in post-swap convention
%                      (row1 = packed-y, row2 = packed-x, row3 = height, row4 = power)
%   antennaObjectMbs : nBands×N antenna template stack
%   extras           : struct array, each entry with fields
%                        .x        world x-coordinate (same frame as generate_hex_sites xs)
%                        .y        world y-coordinate
%                        .height   antenna height (m)
%                        .power    transmit power (W)
%                        .antenna  (optional) 1×nBands antenna template row;
%                                  defaults to defaultAntenna if omitted/empty
%   defaultAntenna   : 1×nBands antenna template row used when an extra omits .antenna
%
% Outputs:
%   mbs_params, antennaObjectMbs : extended matrices
%   containsMbs                  : 1 if any BS present, else 0
%   numMbs                       : total fixed-BS count (base + extras)

    if isempty(extras)
        containsMbs = double(size(mbs_params, 2) > 0);
        numMbs      = size(mbs_params, 2);
        return;
    end

    nBands = size(antennaObjectMbs, 1);
    assert(numel(defaultAntenna) == nBands, ...
        'defaultAntenna must have nBands=%d elements (got %d).', nBands, numel(defaultAntenna));

    nExtra = numel(extras);
    addCols     = zeros(4, nExtra);
    addAntennas = repmat(defaultAntenna(:), 1, nExtra);

    for k = 1:nExtra
        e = extras(k);
        % Match the row1<->row2 swap applied to the base set: pack as
        % [x; y; h; p] then swap rows 1 and 2 so it lines up with mbs_params.
        packed = [e.x; e.y; e.height; e.power];
        tmp = packed(1); packed(1) = packed(2); packed(2) = tmp;
        addCols(:, k) = packed;

        if isfield(e, 'antenna') && ~isempty(e.antenna)
            assert(numel(e.antenna) == nBands, ...
                'extras(%d).antenna must have nBands=%d elements (got %d).', ...
                k, nBands, numel(e.antenna));
            addAntennas(:, k) = e.antenna(:);
        end
    end

    mbs_params       = [mbs_params, addCols];
    antennaObjectMbs = [antennaObjectMbs, addAntennas];
    numMbs           = size(mbs_params, 2);
    containsMbs      = double(numMbs > 0);
end
