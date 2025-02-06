close all;
hold on;
% part a:

N = 100000;
t = linspace(-1.5, 1.5, N);

f_t_exp = exp(-2 * (t.^2));
f_t_sin = sin(5*pi.*t);

f_t = f_t_exp.*f_t_sin; % Pi array
plot(t, f_t);

% part b:
% interval evaluation:

% we will set length(counts) = N / div_factor to arrange the sampled vals
% into sections of the histogram

div_factor = 3000;

bin_centers = linspace(-1, 1, (N / div_factor));

bin_counts = zeros(1, length(bin_centers));
d_f_t = 2/length(bin_centers); % (x_max - x_min) / k


for i = 1:length(f_t)
    val = f_t(i);
    for j = 1:length(bin_centers)
        if( (val > bin_centers(j) - (d_f_t/2)) && (val < bin_centers(j) + (d_f_t/2)) )
            bin_counts(j) = bin_counts(j) + 1;
            break;
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

p_x = (bin_counts./length(bin_counts)./div_factor);
plot(linspace(-1, 1, length(bin_counts)), p_x*100); 

% sum(bin_counts./length(bin_counts)./div_factor) should be around 1 for
% pdf creation 

% part d:

% mean
mean = 0;


for i = 1:length(f_t)
    val = f_t(i);
    ind = min([round(i/div_factor) + 1, length(p_x)]);
    mean = mean + (p_x(ind) * val);
end
mean
% variance

variance = 0;
for i = 1:length(f_t)
    x = f_t(i);
    ind = min([round(i/div_factor) + 1, length(p_x)]);
    ind
    variance = variance + ((x-mean)^2 * p_x(ind));
end

variance


