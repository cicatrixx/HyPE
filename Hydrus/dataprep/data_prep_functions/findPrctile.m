function y = findPrctile(x,p)
    % Return data point nearest to percentile p, without interpolation
    % p:0-100
    % Copied from matlab example:https://nl.mathworks.com/help/matlab/matlab_prog/grouped-calculations-in-tables-and-timetables.html#GroupedCalculationsUsingTimeSeriesDataInTablesExample-10
    xs = sort(x);
    n = sum(~isnan(x)); % use non-NaN elements only
    k = p*n/100 + 0.5;  % index of data point that represents 100*(i-0.5)/n th percentile
    y = xs(round(k));   % data point nearest specified p
end