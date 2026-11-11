function r = competition_rank_rev1(v)
%COMPETITION_RANK_REV1 Rank from the largest, ties sharing the best rank.
%   r = competition_rank_rev1(v) gives 1 to the largest value; tied values
%   share the best rank and the next value skips the tied places
%   ("1224" ranking). Values within a relative 1e-12 count as tied.
%
%   Example:
%     competition_rank_rev1([5 9 9 1])   % [3 1 1 4]
%
%   See also tiedrank.

arguments
    v double {mustBeVector}
end
c = v(:);
r = 1 + sum(c' > c + 1e-12*abs(c), 2);          % 1 + number of values strictly larger
r = reshape(r, size(v));
end
