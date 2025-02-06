close all;
hold on;
% part a:

t = -1.5:0.0001:1.5;
t(end)=[]; % remove end
f_t_exp = exp(-2 * (t.^2));
f_t_sin = sin(5*pi.*t);

f_t = f_t_exp.*f_t_sin; % Pi array
plot(t, f_t);


% part b:
% interval evaluation:

% we will set length(counts) = N / div_factor to arrange the sampled vals
% into sections of the histogram
num_cols = 300;
div_factor = length(t)/num_cols;

bin_centers = linspace(-1, 1, num_cols);

bin_counts = zeros(1, length(bin_centers));

d_f_t = 2/length(bin_centers); % (x_max - x_min) / k

% n = bin count in a bin (bin_counts(j))

% N = sample size (length(f_t))

bin_centers = bin_centers+(d_f_t/2);

for i = 1:length(f_t)
    val = f_t(i);
    for j = 1:length(bin_centers)
        if( (val > bin_centers(j) - (d_f_t/2)) && (val < bin_centers(j) + (d_f_t/2)) )
            bin_counts(j) = bin_counts(j) + 1;
            
        end
    end
end

figure;
bar(linspace(-1, 1, length(bin_counts)), bin_counts)

% part c:

% now what we want to plot is the relative probability of the values.
% to do this we can simply normalize the bin center counts by the number of
% bins
figure;


p_x = bin_counts / sum(bin_counts);
plot(linspace(-1, 1, length(p_x)), p_x*100); 

sum(p_x);
% pdf creation 

% part d:

% mean
mean = sum(bin_centers .* p_x * d_f_t)
% variance

variance = 0;
test = [];
for i = 1:length(f_t)
    x = f_t(i);
    ind = min([floor(i/div_factor) + 1, length(p_x)]);
    test(ind) = p_x(ind);
    variance = variance + ((x-mean)^2 * p_x(ind));
end

variance;


