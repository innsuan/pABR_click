function data_out = draw_pABR_clicks()
% YW @ 05/06/2026
% Change path_read and name_read to match the example data location
% Preprocessed data include y, fs, and list
path_read = '';% set to folder containing the data file
name_read = 'B323011003.mat';
load([path_read name_read]);
% Loaded variables:
%   y    - neural recording, [n_epochs x pts_epoch] in uV
%          each row is one epoch (one full Poisson pip train presentation)
%   fs   - sampling frequency (Hz)
%   list - stimulus parameter structure with fields:
%            list.level   - [n_epochs x 1] dB SPL for each epoch
%            list.i_onset - [n_epochs x n_clicks_per_epoch] pip onset indices
%                           within each epoch (in samples).
%                           IMPORTANT: sign of i_onset encodes polarity —
%                           positive index = condensation click (+1)
%                           negative index = rarefaction click (-1)
%                           use abs() for indexing, sign() for polarity

all_levels = unique(list.level);
n_levels = length(all_levels);
n_epochs = size(y,1);
all_polarities = [-1 1];
pts_raw = round(100e-3*fs);% 100 ms
pts_abr = round(30e-3*fs);% 30 ms
t = (0:(pts_raw - 1))/fs;% ms
t = t*1e3;%
% list of parameters for each individual click
n_pips_per_epoch = size(list.i_onset, 2);
n_clicks = n_epochs * n_pips_per_epoch;
%% parsing data
list_ici = zeros(1, n_clicks);
list_polarity_2 = zeros(1, n_clicks);
list_level = zeros(1, n_clicks);
y_click = zeros(n_clicks, pts_raw);

for i_epoch = 1:n_epochs
    i_start = (i_epoch-1)*n_pips_per_epoch + 1;
    i_end = i_epoch*n_pips_per_epoch;
    i_onset_this = abs(list.i_onset(i_epoch,:));
    list_ici(i_start:i_end) = [Inf diff(i_onset_this)];
    list_polarity_2(i_start:i_end) = sign(list.i_onset(i_epoch,:));
    list_level(i_start:i_end) = list.level(i_epoch);
    for ii_onset = 1:n_pips_per_epoch
        y_click(i_start+ii_onset-1,:) = y(i_epoch, i_onset_this(ii_onset):i_onset_this(ii_onset)+pts_raw-1);
    end
end
list_polarity_1 = [nan list_polarity_2];
list_polarity_1 = list_polarity_1(1:length(list_polarity_2));
list_ici = list_ici/fs*1000;% in ms

%% get pABR grand average
abr_wav = zeros(n_levels, pts_raw);
cm_wav = zeros(n_levels, pts_raw);
abr_wav_raw = cell(n_levels, 1);
for i_level = 1:n_levels
    i_this = list_level == all_levels(i_level);  % was ispl_list == i_level
    abr_wav_raw{i_level} = y_click(i_this, 1:pts_raw);
    abr_wav(i_level,:) = mean(abr_wav_raw{i_level}, 'omitnan');
    for i_pol = 1:length(all_polarities)
        i_this_pol = i_this & list_polarity_2 == all_polarities(i_pol);  % was ipol_list == pols(i_pol)
        abr_wav_pol(i_level, i_pol,:) = mean(y_click(i_this_pol, 1:pts_raw), 'omitnan');
    end
    cm_wav(i_level,:) = (squeeze(abr_wav_pol(i_level,1,:)) - squeeze(abr_wav_pol(i_level,2,:)))'/2;
end

%% figure for pABR grand average
fig_grand = figure('Position',[0 0 500 900],...
    'CreateFcn',{@movegui,'northwest'});
figure(fig_grand);
n_rows = n_levels;
ylims = zeros(n_levels, 2);
count = 0;
for i_level = n_levels:-1:1
    count = count + 1;
    h_grand(i_level) = subplot(n_rows,1,count);
    hold on
    plot(t(1:pts_abr),squeeze(abr_wav(i_level,1:pts_abr)),'k');
    plot(t(1:pts_abr),squeeze(cm_wav(i_level,1:pts_abr)),'b');
    ylims(count,:) = get(gca,'YLim');
    grid on
end
ylims = [min(ylims(:,1)) max(ylims(:,2))];
set(h_grand,'YLim',ylims)
for i_level = 1:n_levels
    text(h_grand(i_level),mean(t(1:pts_abr)),ylims(2),[num2str(all_levels(i_level)) 'dB'],...
        'HorizontalAlignment','center','VerticalAlignment','bottom')
end

%% ICI analysis
abr_wav_diff = [];
lambda = mean(list_ici(list_ici<350));
for i_level = 1:n_levels
    ici_window_init = 1; % in ms
    ici_min = 0.5;
    % ici_step = 1;
    ici_max = 75;%lambda*3;
    ici_nsteps = 75;
    i = 0;
    abr_wav_prev_n = squeeze(abr_wav_pol(i_level,1,:))';
    abr_wav_prev_p = squeeze(abr_wav_pol(i_level,2,:))';
    all_icis = logspace(log10(ici_min),log10(ici_max),ici_nsteps);%[ici_min:0.25:10 logspace(log10(10),log10(ici_max),ici_nsteps-length(ici_min:0.25:10))];% logspace(log10(ici_min),log10(ici_max),ici_nsteps);%ici_min:ici_step:ici_max;%[1 2 5 10 20 50 100];%logspace(log10(ici_min),log10(ici_max),ici_nsteps);%ici_min:ici_step:ici_max;
    all_icis = unique(all_icis);
    abr_wav_diff{i_level} = zeros(length(all_icis),pts_raw);
    % first iteration
    ici_this_level = list_ici(list_level == all_levels(i_level));
    polarity_2_this_level = list_polarity_2(list_level == all_levels(i_level));
    polarity_1_this_level = list_polarity_1(list_level == all_levels(i_level));
    for current_ici = all_icis
        i = i + 1;
        if current_ici == ici_min
            ici_window = ici_window_init;
            window_param = sqrt(ici_window^2/(current_ici^2*4) + 1) + ici_window/2/current_ici;
            p_in_window = expcdf(current_ici*window_param,lambda) - expcdf(current_ici/window_param,lambda);%
        else
            window_param = fzero(@(x) expcdf(current_ici*x,lambda) - expcdf(current_ici/x,lambda) - p_in_window,ici_window);
        end
        ici_lower(i) = current_ici/window_param;
        ici_upper(i) = current_ici*window_param;
        i_this = ici_this_level > ici_lower(i) & ici_this_level < ici_upper(i) ...
            & ~(ici_this_level == Inf);
        n_ici(i) = sum(i_this);
        abr_wav_this = abr_wav_raw{i_level}(i_this,:);
        ici_this = ici_this_level(i_this);
        polarity_2_this = polarity_2_this_level(i_this);
        polarity_1_this = polarity_1_this_level(i_this);
        abr_wav_diff_raw = nan(length(ici_this),pts_raw);
        for ii = 1:length(ici_this)
            if polarity_1_this(ii) == -1
                abr_wav_prev = abr_wav_prev_n;
            elseif polarity_1_this(ii) == 1
                abr_wav_prev = abr_wav_prev_p;
            else
                abr_wav_prev = (abr_wav_prev_n + abr_wav_prev_p)/2;
            end
            ici_pts = round(ici_this(ii)/1000*fs);
            abr_wav_diff_this = [zeros(1,ici_pts) abr_wav_this(ii,:)]...
                - [abr_wav_prev zeros(1,ici_pts)];
            abr_wav_diff_raw(ii,:) = abr_wav_diff_this(ici_pts+1:ici_pts+pts_raw);
        end
        % least square linear fit
        x_lm = ici_this;
        y_lm = abr_wav_diff_raw;
        cm_wav_diff{i_level}(i,:) = (mean(abr_wav_diff_raw(polarity_2_this == -1,:),'omitnan') - mean(abr_wav_diff_raw(polarity_2_this == 1,:),'omitnan'))/2;
        for i_sample = 1:size(y_lm,2)
            p = polyfit(x_lm((~isnan(y_lm(:,i_sample))))',y_lm((~isnan(y_lm(:,i_sample))),i_sample),1);
            c_1(i_sample) = p(1);
            c_2(i_sample) = p(2);
        end
        abr_wav_diff{i_level}(i,:) = c_2 + c_1*current_ici; %mean(abr_wav_diff_raw,'omitnan') - mean(mean(abr_wav_diff_raw(:,end-pts_abr:end),'omitnan'));
    end
    icis{i_level} = all_icis;%mean(ici_this);
end

%% figure for ici plotting
f_ici = figure('unit','inch','Position',[0 0 3 9]);
for i_level = 1:n_levels
    h_ici(i_level) = subplot(1,n_levels,i_level);
end
gap_sample = 3;
gap_v = max(abs(abr_wav_diff{end}(:)))/3/gap_sample;
for i_level = 1:n_levels
    y_baseline = abr_wav_diff{i_level}(:,round(50/1000*fs):round(50/1000*fs)+pts_abr-1);
    y_abr = abr_wav_diff{i_level};
    y_cm = cm_wav_diff{i_level};
    axes(h_ici(i_level));hold on;
    for i_ici = 1:gap_sample:length(all_icis)
        plot(t(1:pts_abr),y_baseline(i_ici,1:pts_abr)+i_ici*gap_v,'color',[0.9 0.9 0.9]);
        plot(t(1:pts_abr),y_abr(i_ici,1:pts_abr)+i_ici*gap_v,'k');
        % plot(t(1:pts_abr),y_cm(i_ici,1:pts_abr)+i_ici*gap_v,'b');
    end
    set(gca,'Ytick',[],'YTickLabel',[],...
        'xtick',[2 5 7 10 12 15],'xTickLabel',[2 5 7 10 12 15]);
    xlim([0 20]);
    ylim([-gap_v (i+3*gap_sample)*gap_v]);
    grid on;
    % end
end

%% data output
data_out.t = t;
data_out.icis = icis;
data_out.all_levels = all_levels;
data_out.y_abr = abr_wav;
data_out.y_ici = abr_wav_diff;
end