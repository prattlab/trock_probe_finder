%% Preamble
% This script aims to develop a matched filter for detection of optimal
% drop-off and reference probes along TP53 based on data describing
% frequency of mutations for different patients. The script will be
% rewritten in R but is starting in MATLAB to allow a quick demo.

% Edwin Rock, BME PhD Candidate rotating in Erica Pratt Lab at Boston
% University. 9/29/24

%% Parameters
probe_min = 19; % Shortest length each probe can be
probe_max = 31; % Longest length each probe can be
space_max = 30; % Longest space in between probes
ref_probe_size = 21; % Allows reference probes to be pre-allocated at a certain size
a = 3; % Parameter affecting score calculation

%% Data Setup

% Reading table from CSV and isolating vector of amino acid locations with
% mutations associated
T = readtable("cBioPortal_data.csv", 'NumHeaderLines', 1);
mutation_locations = T{:, 3};

% Multiplying locations by 3 to change from amino acids to nucleotide bases
mutation_locations = mutation_locations*3;

% Converting mutation_locations into a vector mutations with number of
% mutations at each location (find a nicer way to do this in the future)
mutations = zeros(max(mutation_locations), 1);
for i = 1:max(mutation_locations)
    mutation_present = mutation_locations == i;
    mutations(i) = sum(mutation_present);
end

%% Matched Filter Reward Vector

% Creating a vector with high value for the middle points (minimum
% probe length) and lower value for distance away from the middle (moving
% towards maximum probe length
ends = 1:(probe_max-probe_min);
reward = [ends (max(ends)+1)*ones(1, probe_min) flip(ends)]';

%% Convolution

% Convolving reward vector with mutations vector to see highest/lowest
% association
association = conv(mutations, reward);

% Creating vector conv_location which matches each point in association to
% the location of the center base pair of the probe
conv_location = (2-ceil(length(reward)/2)):(length(mutations)+floor(length(reward)/2));

%% Data Visualization

figure(1)

% Lollipop plot of mutation locations
subplot(2, 2, 1)
mutations_plot = mutations;
mutations_plot(mutations == 0) = NaN;
stem(mutations_plot, ...
    'LineWidth', 1, ...
    'Color', 'k', ...
    'Marker', 'none')
axis('padded')
xlabel('Nucleotide Base', ...
    'FontSize', 18)
ylabel('# of Mutations', ...
    'FontSize', 18)
title('Lollipop Plot of Mutations along TP53', ...
    'FontSize', 20)

% Reward vector
subplot(2, 2, 2)
stairs(0:(length(reward)+2), [0; reward; 0; 0], ...
    'LineWidth', 2, ...
    'Color', 'k')
xlim([0 length(reward)+2])
ylim([0 max(reward)+1])
xlabel('Probe Length', ...
    'FontSize', 18)
ylabel('Reward', ...
    'FontSize', 18)
title('Reward Vector', ...
    'FontSize', 20)
grid on

% Output
subplot(2, 1, 2)
plot(conv_location, association, ...
    'LineWidth', 2, ...
    'Color', 'k')
xlim([conv_location(1) conv_location(end)])
ylim([0 max(association)*1.2])
xlabel('Probe Center', ...
    'FontSize', 18)
ylabel('Reward', ...
    'FontSize', 18)
title('Association of Lollipop Plot with Reward', ...
    'FontSize', 20)
xline(1)
xline(length(mutations))
grid on

% Figure Details
set(gcf, 'Color', 'w')

%% Identifying Optimized Probe Center Pairs

% Finding maximum value along association vector
% M is value, I is index at which it occurs
[~, I_DO] = max(association);

% Finding window around drop-off probe to look for center of reference
% probe and checking whether this window goes out of the length of the gene
open_window = probe_min - 1 + 2*(probe_max - probe_min) + space_max;
lower_bound = I_DO - open_window;
upper_bound = I_DO + open_window;
if lower_bound <= 0
    lower_bound = 1;
end
if upper_bound > length(association)
    upper_bound = length(association);
end
[~, I_R] = min(association(lower_bound:upper_bound));

% Converting I to location of probe center
I_DO = conv_location(I_DO);
I_R = conv_location(I_R + lower_bound - 1);

%% Finding Ideal Probe Length
% NOTE: FIND A MORE EFFICIENT WAY TO DO THIS

% probe_extensions holds each probe's optimal extension on the left and
% right in a row
probe_extensions = zeros(length(I_DO_R), 2);

% For each index in I_DO, extend along either side and see how this changes
% the score to maximize score
for i = 1:length(I_DO_R)
    
    % left_right is a matrix with number of extra bases on left as columns
    % [0; 1; 2; 3 ... n] and extra bases on right as rows [0 1 2 3 ... n].
    % Values are the score for that setup
    left_right = zeros(probe_max - probe_min + 1);

    % j is bases on left side, k is on right side
    for j = 0:(probe_max-probe_min)
        for k = 0:(probe_max-probe_min)
            
            % Checking still under probe_max
            if j + k > (probe_max-probe_min)
                continue
            end

            % Scoring configuration
            left_right(j+1, k+1) = sum(mutations((I_DO_R(i)-floor(probe_min/2)-j):(I_DO_R(i)+floor(probe_min/2)+k)))/(a*sum(mutations((I_R_R(i)-floor(probe_min/2)):(I_R_R(i)-floor(probe_min/2))))+(probe_min+j+k)*(ref_probe_size)*abs(j+floor(probe_min/2)+abs(I_DO_R(i)-I_R_R(i))+1+floor(ref_probe_size)-89.9));
        end
    end

    % Saving optimal values to probe_extensions
    [x,y] = find(left_right == max(left_right, [], 'all'));
    probe_extensions(i, :) = [x-1 y-1];
    
end

%% Visualizing Probes Along Lollipop Plot

I_DO_R = [final_probes.drop_off_start(1:7) final_probes.drop_off_end(1:7)];
I_R_R = [final_probes.ref_start(1:7) final_probes.ref_end(1:7)];

figure(3)
stem(mutations_plot, ...
    'LineWidth', 1, ...
    'Color', 'k', ...
    'Marker', 'none')
hold on
for i = 1:length(I_DO_R)
    stem((I_DO_R(i, 1)):(I_DO_R(i, 2)), mutations_plot((I_DO_R(i, 1)):(I_DO_R(i, 2))), ...
        'LineWidth', 1, ...
        'Color', 'b', ...
        'Marker', 'none')
    stem((I_R_R(i, 1)):(I_R_R(i, 2)), mutations_plot((I_R_R(i, 1)):(I_R_R(i, 2))), ...
        'LineWidth', 1, ...
        'Color', 'r', ...
        'Marker', 'none')
end
hold off
axis('padded')
xlabel('Nucleotide Base', ...
    'FontSize', 18)
ylabel('# of Mutations', ...
    'FontSize', 18)
title('Lollipop Plot of Mutations along TP53', ...
    'FontSize', 20)

%% Another visualization
figure(4)

% subplot(2, 1, 1)
% stem(mutations_plot, ...
%     'LineWidth', 1.5, ...
%     'Color', 'k', ...
%     'Marker', 'none')
% hold on
% for i = 1:length(I_DO_R)
%     rectangle(Position=[I_DO_brute(i, 1) 0 (I_DO_brute(i, 2) - I_DO_brute(i, 1) + 1) 1.2*max(mutations_plot)], ...
%         FaceColor = [0 0 1 0.3], ...
%         EdgeColor = [0 0 1 0.3])
%     rectangle(Position=[I_R_brute(i, 1) 0 (I_R_brute(i, 2) - I_R_brute(i, 1) + 1) 1.2*max(mutations_plot)], ...
%         FaceColor = [1 0 0 0.3], ...
%         EdgeColor = [1 0 0 0.3])
% end
% hold off
% ylim([0 1.1*max(mutations_plot)])
% title('Brute Force Approach', ...
%     'FontSize', 20)
% xlabel('Nucleotide Base', ...
%     'FontSize', 18)
% ylabel('# of Mutations', ...
%     'FontSize', 18)
% 
% subplot(2, 1, 2)
subplot(2, 1, 1)
stem(mutations_plot, ...
    'LineWidth', 1.5, ...
    'Color', 'k', ...
    'Marker', 'none')
hold on
for i = 1:length(I_DO_R)
    rectangle(Position=[I_DO_R(i, 1) 0 (I_DO_R(i, 2) - I_DO_R(i, 1) + 1) 1.2*max(mutations_plot)], ...
        FaceColor = [0 0 1 0.3], ...
        EdgeColor = [0 0 1 0.3])
    rectangle(Position=[I_R_R(i, 1) 0 (I_R_R(i, 2) - I_R_R(i, 1) + 1) 1.2*max(mutations_plot)], ...
        FaceColor = [1 0 0 0.3], ...
        EdgeColor = [1 0 0 0.3])
end
hold off
ylim([0 1.1*max(mutations_plot)])
title('Probes Across Entire Gene', ...
    'FontSize', 20)
ylabel('# of Mutations', ...
    'FontSize', 18)
grid on
d_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color',[0 0 1 0.3]);
r_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color', [1 0 0 0.3]);
legend([d_line r_line], ...
    {'Drop Off Probes', 'Reference Probes'}, ...
    'FontSize', 12)

subplot(2, 1, 2)
stem(mutations_plot, ...
    'LineWidth', 1.5, ...
    'Color', 'k', ...
    'Marker', 'none')
hold on
for i = 1:length(I_DO_R)
    rectangle(Position=[I_DO_R(i, 1) 0 (I_DO_R(i, 2) - I_DO_R(i, 1) + 1) 1.2*max(mutations_plot)], ...
        FaceColor = [0 0 1 0.3], ...
        EdgeColor = [0 0 1 0.3])
    rectangle(Position=[I_R_R(i, 1) 0 (I_R_R(i, 2) - I_R_R(i, 1) + 1) 1.2*max(mutations_plot)], ...
        FaceColor = [1 0 0 0.3], ...
        EdgeColor = [1 0 0 0.3])
end
hold off
xlim([4000 5800])
ylim([0 1.1*max(mutations_plot)])
title('Probes Across Mutational Hotspots', ...
    'FontSize', 20)
xlabel('Nucleotide Base', ...
    'FontSize', 18)
ylabel('# of Mutations', ...
    'FontSize', 18)
grid on
d_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color',[0 0 1 0.3]);
r_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color', [1 0 0 0.3]);
legend([d_line r_line], ...
    {'Drop Off Probes', 'Reference Probes'}, ...
    'FontSize', 12)

set(gcf, 'Color', 'w')

%%

figure(4)
hold on
for i = 1:length(I_DO_R)
    rectangle(Position=[I_DO_R(i)-10, 0, 21,1.3*max(mutations)], ...
        FaceColor='#0072BD', ...
        EdgeColor='k')
    rectangle(Position=[I_R_R(i)-10, 0, 21,1.3*max(mutations)], ...
        FaceColor='#D95319', ...
        EdgeColor='k')
end
stem(mutations_plot, ...
    'LineWidth', 1, ...
    'Color', 'k', ...
    'Marker', 'none')
hold off
%xlim([-100 1300])
ylim([0 1.3*max(mutations)])
xlabel('Nucleotide Base', ...
    'FontSize', 18)
ylabel('# of Mutations', ...
    'FontSize', 18)
title('Lollipop Plot of Mutations along TP53', ...
    'FontSize', 20)

%%

triangle_pulse = [zeros(1,2) (0:6)/3 zeros(1,3)];
square_pulse = [zeros(1,14) 2*ones(1,5) zeros(1,10)];

figure(5)

subplot(3, 2, 1)
stairs(square_pulse, ...
    'Color', '#0072BD', ...
    'LineWidth', 2)
grid on
axis padded
title('Square Pulse', ...
    'FontSize', 13)

subplot(3, 2, 2)
plot(triangle_pulse, ...
    'Color', '#D95319', ...
    'LineWidth', 2)
grid on
axis padded
title('Triangle Pulse', ...
    'FontSize', 13)

subplot(3, 1, 2)
plot(conv(square_pulse, triangle_pulse), ...
    'color', '#7E2F8E', ...
    'LineWidth', 2)
grid on
axis padded
title('Convolution Output', ...
    'FontSize', 13)

subplot(3, 3, 7)
stairs(square_pulse, ...
    'color', '#0072BD', ...
    'LineWidth', 2)
hold on
plot(-5:6, flip(triangle_pulse), ...
    'Color', '#D95319', ...
    'LineWidth', 2)
hold off
grid on
axis padded

subplot(3, 3, 8)
stairs(square_pulse, ...
    'color', '#0072BD', ...
    'LineWidth', 2)
hold on
plot(10:21, flip(triangle_pulse), ...
    'Color', '#D95319', ...
    'LineWidth', 2)
hold off
grid on
axis padded

subplot(3, 3, 9)
stairs(square_pulse, ...
    'color', '#0072BD', ...
    'LineWidth', 2)
hold on
plot(23:34, flip(triangle_pulse), ...
    'Color', '#D95319', ...
    'LineWidth', 2)
hold off
grid on
axis padded

set(gcf, 'Color', 'w')

% Matched filter demo

t = 1:180;
%signal_alone = awgn(ones(size(t)), 1e-4);
signal_alone = [zeros(1, 10) ones(1, 25) -1.^(1:25) -1*ones(1, 25) 0 sin(0.25*(1:38))];
signal_pulse = awgn([zeros(1, 800) signal_alone zeros(1, 600)], -1);

figure(6)

subplot(2, 2, 1)
plot(signal_pulse, ...
    'Color', '#0072BD', ...
    'LineWidth', 2)
grid on
axis padded
title('Noisy Signal', 'FontSize', 18)

subplot(2, 2, 2)
plot(signal_alone, ...
    'Color', '#D95319', ...
    'LineWidth', 2)
grid on
axis padded
title('Our Signal', 'FontSize', 18)

subplot(2, 1, 2)
plot(conv(signal_pulse, flip(signal_alone)), ...
    'Color', '#7E2F8E', ...
    'LineWidth', 2)
grid on
axis padded
title('Matched Filter Output', 'FontSize', 18)

set(gcf, 'Color', 'w')

%% Plotting Two Different Probe Sets

figure;

yline(1); yline(2)
ylim([0 3])
xlim([4000 7000])

for i = 1:7
    rectangle(Position=[final_probes_train.drop_off_start(i) 0.6 (final_probes_train.drop_off_end(i)-final_probes_train.drop_off_start(i)) 0.8], ...
            FaceColor = [0 0 1 1], ...
            EdgeColor = [0 0 1 1])
    rectangle(Position=[final_probes_train.ref_start(i) 0.6 (final_probes_train.ref_end(i)-final_probes_train.ref_start(i)) 0.8], ...
        FaceColor = [1 0 0 1], ...
        EdgeColor = [1 0 0 1])
    rectangle(Position=[final_probes_test.drop_off_start(i) 1.6 (final_probes_test.drop_off_end(i)-final_probes_test.drop_off_start(i)) 0.8], ...
        FaceColor = [0 0 1 1], ...
        EdgeColor = [0 0 1 1])
    rectangle(Position=[final_probes_test.ref_start(i) 1.6 (final_probes_test.ref_end(i)-final_probes_test.ref_start(i)) 0.8], ...
        FaceColor = [1 0 0 1], ...
        EdgeColor = [1 0 0 1])
end

set(gcf, 'Color', 'w')
set(gca, 'ytick', [])

ylabel('Train vs. Test', 'FontSize', 14)
xlabel('Nucleotide Base', 'FontSize', 14)
title('Probe Locations from Split Datasets', 'FontSize', 14)
set(gca, 'box', 'on')


%%

figure;

mut_1 = train_mutations_plot;
mut_2 = test_mutations_plot;

assay_1 = [final_probes_train.drop_off_start final_probes_train.drop_off_end final_probes_train.ref_start final_probes_train.ref_end];
%assay_2 = [final_probes_test.drop_off_start final_probes_test.drop_off_end final_probes_test.ref_start final_probes_test.ref_end];
assay_2 = assay_1;

subplot(2, 1, 1)
stem(mut_1, ...
    'LineWidth', 1.5, ...
    'Color', 'k', ...
    'Marker', 'none')
hold on
for i = 1:size(assay_1, 1)
    rectangle(Position=[assay_1(i, 1) 0 (assay_1(i, 2) - assay_1(i, 1) + 1) 1.2*max(mut_1)], ...
        FaceColor = [0 0 1 0.3], ...
        EdgeColor = [0 0 1 0.3])
    rectangle(Position=[assay_1(i, 3) 0 (assay_1(i, 4) - assay_1(i, 3) + 1) 1.2*max(mut_1)], ...
        FaceColor = [1 0 0 0.3], ...
        EdgeColor = [1 0 0 0.3])
end
hold off
xlim([3800 6800])
ylim([0 1.1*max(mut_1)])
title('Probes Across Train Dataset', ...
    'FontSize', 20)
ylabel('# of Mutations', ...
    'FontSize', 18)
grid on
d_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color',[0 0 1 0.3]);
r_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color', [1 0 0 0.3]);
legend([d_line r_line], ...
    {'Drop Off Probes', 'Reference Probes'}, ...
    'FontSize', 14)

subplot(2, 1, 2)
stem(mut_2, ...
    'LineWidth', 1.5, ...
    'Color', 'k', ...
    'Marker', 'none')
hold on
for i = 1:length(I_DO_R)
    rectangle(Position=[assay_2(i, 1) 0 (assay_2(i, 2) - assay_2(i, 1) + 1) 1.2*max(mut_2)], ...
        FaceColor = [0 0 1 0.3], ...
        EdgeColor = [0 0 1 0.3])
    rectangle(Position=[assay_2(i, 3) 0 (assay_2(i, 4) - assay_2(i, 3) + 1) 1.2*max(mut_2)], ...
        FaceColor = [1 0 0 0.3], ...
        EdgeColor = [1 0 0 0.3])
end
hold off
xlim([3800 6800])
ylim([0 1.1*max(mut_2)])
title('Probes Across Test Dataset', ...
    'FontSize', 20)
xlabel('Nucleotide Base', ...
    'FontSize', 18)
ylabel('# of Mutations', ...
    'FontSize', 18)
grid on
d_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color',[0 0 1 0.3]);
r_line = line(NaN, NaN, ...
    'LineWidth', 12, ...
    'Color', [1 0 0 0.3]);
legend([d_line r_line], ...
    {'Drop Off Probes', 'Reference Probes'}, ...
    'FontSize', 14)

set(gcf, 'Color', 'w')


%%

triangle_pulse = [zeros(1,2) (0:6)/3 zeros(1,3)];
square_pulse = [zeros(1,14) 2*ones(1,5) zeros(1,10)];

conv_length = length(triangle_pulse) + length(square_pulse) - 1;
conv_out = conv(square_pulse, triangle_pulse);

writerObj = VideoWriter('Convolution', 'MPEG-4');
writerObj.FrameRate = 7.5;
open(writerObj)

figure;

for i = 1:conv_length
    subplot(2, 1, 1)
    hold on
    stairs(1:length(square_pulse), square_pulse, ...
        'Color', '#0072BD', ...
        'LineWidth', 2)
    triangle = plot(flip(i+1-(1:length(triangle_pulse))), flip(triangle_pulse), ...
        'Color', '#D95319', ...
        'LineWidth', 2);
    hold off
    grid on
    title('Signal Lineup', ...
        'FontSize', 13)
    xlim([-11 41])

    subplot(2, 1, 2)
    plot(1:i, conv_out(1:i), ...
        'Color', '#7E2F8E', ...
        'LineWidth', 2)
    grid on
    title('Convolution Output', ...
        'FontSize', 13)
    xlim([1 40])
    
    %pause(0.25)
    set(gcf, 'Color', 'w')
    frame = getframe(gcf);
    writeVideo(writerObj, frame)
    
    if i < conv_length
        delete(triangle)
    end
end

close(writerObj)


