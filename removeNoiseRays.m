function good_indices = removeNoiseRays(data,x,y,z, thre_snr, varargin)
% For Siemens reconstruction
% Remove the rays with noise spike, and update the trajectory matrix
if nargin == 5
    tail = 10;
elseif nargin == 6
    tail = varargin{1};
else
    ME = MException('not a valid number of arguments');
    throw(ME);
end
thre = thre_snr*mean(mean(abs(data(1:1:5,:))));
data_tail = abs(data(tail+1:end,:));
max_tail = max(data_tail);
good_indices = find(max_tail < thre);

end

