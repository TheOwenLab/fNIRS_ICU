function DataTemp = DiscardDataPoints(DataTemp, varargin)

if nargin < 2
    discard_flag = 1;
    h = figure;
    while discard_flag
        [value position] = max(mean(DataTemp.d,1));
        plot(DataTemp.s), hold on
        plot(DataTemp.d(:,position(1))./max(DataTemp.d(:,position(1))));
        hold off
        discard_range = input('Any data points to discard (e.g., 1:154)? Enter to continue: ');
        if ~isempty(discard_range)
            DataTemp.t(discard_range) = [];
            DataTemp.d(discard_range,:) = [];
            DataTemp.s(discard_range,:) = [];
            DataTemp.aux(discard_range,:) = [];
        else
            discard_flag = 0;
        end
    end
else
    DataTemp.t(varargin{1}) = [];
    DataTemp.d(varargin{1},:) = [];
    DataTemp.s(varargin{1},:) = [];
    DataTemp.aux(varargin{1},:) = [];
    if ~isempty(DataTemp.dc)
        DataTemp.dc(varargin{1},:,:) = [];
    end
end

return