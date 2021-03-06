%
% Mixer waveform generator
% Sweep [SIGNAL] from noise_freq_start ~ noise_freq_end
%

function y = mixer(desc, test_mode, noise_freq, signal_freq_start, signal_freq_end, signal_freq_step, signal_phase_cnt, noise_level, int_threshold)
%function y = mixer(test_mode, mixer_freq, noise_freq_start, noise_freq_end, noise_freq_step, noise_phase_cnt, noise_level, int_threshold)

% clear
close all
% clc

fs = 16e6;   % 16M
mixing_duration = 200e-6;   % 200us

noise_type = ['Sina    '; 'Sqr     '; 'Sina*H  '; 'Sqr*H   '; 'Sina*H*H'; 'Sqr*H*H '];
plot_color = ['b'; 'r'; 'c'; 'k'; 'g'; 'm'];

noise_types = size(noise_type);
noise_types = noise_types(1);
            % BIT0: Sina/Sqr Signal with Sina Noise
            % BIT1: Sina/Sqr Signal with Sqr  Noise
            % BIT2: Sina/Sqr Signal with Sina Noise * Hann
            % BIT3: Sina/Sqr Signal with Sqr  Noise * Hann
            % BIT4: Sina/Sqr Signal with Sina Noise * Hann * Hann
            % BIT5: Sina/Sqr Signal with Sqr  Noise * Hann * Hann
measure_types = 4;
            % MAX_INT
            % MAX_INT_PHASE
            % MIN_INT
            % MIN_INT_PHASE

noise_type_str = cellstr(noise_type);
save_ws = 1

if nargin ~= 9
    disp('[Usage]:')
    disp('        Input1: test_mode')
    disp('        Input2: noise_freq')
    disp('        Input3: sinal_freq_start')
    disp('        Input4: sinal_freq_end')
    disp('        Input5: sinal_freq_step')
    disp('        Input6: sinal_phase_cnt')
    disp('        Input7: noise_level')
    disp('        Input8: threshold')
    % Default setting
    test_mode         = (2^noise_types) - 1;    % Input 2, Select all bits
    noise_freq        = 32.75e3;                % Input 3, 32.75K
    signal_freq_start = 1e3;                    % Input 4, 1K
    signal_freq_end   = 2e6;                    % Input 5, 2M
    signal_freq_step  = 200;                    % Input 6, 0.2K
    signal_phase_cnt  = 36;                     % Input 7, 360/36 = 10
    noise_level       = 1;                      % Input 8
    int_threshold     = 0;                      % Input 9
    if nargin == 1
        % Run on default fully setting
        desc          = sprintf('full_%s', desc);
    else
        % No save
        desc          = 'NoSave';
        save_ws       = 0;
    end
end

%Check Input Parameters
if test_mode == 0 || test_mode >= (2^noise_types)
    disp('Bad test_mode');
    return;
end

signal_freq_cnt     = (signal_freq_end - signal_freq_start) / signal_freq_step;   % (1000 - 40)/10
signal_phase_offset = (2 * pi) / signal_phase_cnt;

disp(        '....................................................................');
disp(        '............................ Test case .............................');
disp(sprintf('Noise freq %d Khz', noise_freq/1000));
disp(sprintf('Sweep noise from %d Khz ~ %d Khz [step: %d Hz] [phase_cnt: %d]', signal_freq_start/1000, signal_freq_end/1000, signal_freq_step, signal_phase_cnt));
disp(sprintf('      noise level %d', noise_level));
disp(sprintf('      check INT Threshold %d', int_threshold));
disp(sprintf('      desc: %s', desc));
test_list = [];
for curr_noise=1:(noise_types)
    if bitand(test_mode, 2^(curr_noise-1))
        test_list = [test_list; noise_type_str(curr_noise)];
    end
end
disp(test_list)
disp(        '....................................................................');


NFFT         = 2^17;
n_samples    = round(mixing_duration*fs);

% Result
test_freq         = zeros(1,signal_freq_cnt);
curr_log_sina     = zeros(noise_types, measure_types);    % element-1: max_data, ele-2: max_phase, ele-3: min, ele-4: min_phase
curr_log_sqr      = zeros(noise_types, measure_types);    % element-1: max_data, ele-2: max_phase, ele-3: min, ele-4: min_phase
sina_max_log      = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sina_maxphase_log = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sina_min_log      = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sina_minphase_log = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sqr_max_log       = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sqr_maxphase_log  = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sqr_min_log       = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4
sqr_minphase_log  = zeros(noise_types, signal_freq_cnt);   % element-1: noise-1, ele-2: noise-2, ele-3: noise-3, ele-4: noise-4


t             = 0:1/fs:mixing_duration-1/fs;
sina_noise    = sin(2*pi*noise_freq*t) * noise_level;
sqr_noise     = square(2*pi*noise_freq*t) * noise_level;
hann_win      = hann(n_samples)';
f             = linspace(0,fs/2,NFFT/2);
% plot(t, hann_win, 'b');
hann_hann_win = hann_win .* hann_win;


for signal_step=0:signal_freq_cnt
% for noise_step=0:2
    signal_freq = signal_freq_start + signal_freq_step * signal_step;

    if rem(signal_freq, 100e3) == 0
        disp(sprintf('Sweep to %d Khz', signal_freq/1000));
    end

    %
    % Generate noise from noise_freq_start to noise_freq_end
    %
    curr_log_sina = zeros(noise_types, measure_types);
    curr_log_sqr  = zeros(noise_types, measure_types);

    for phase_step=0:(signal_phase_cnt-1)
        %
        % Sina Signal
        %
        sina_signal   = sin(2*pi*signal_freq*t + signal_phase_offset * phase_step);
        %
        % Sqr Signal
        %
        square_signal = square(2*pi*signal_freq*t + signal_phase_offset * phase_step);
        for curr_noise=1:(noise_types)
            % BIT0: Sina/Sqr Signal with Sina Noise
            % BIT1: Sina/Sqr Signal with Sqr  Noise
            % BIT2: Sina/Sqr Signal with Sina Noise * Hann
            % BIT3: Sina/Sqr Signal with Sqr  Noise * Hann
            % BIT4: Sina/Sqr Signal with Sina Noise * Hann * Hann
            % BIT5: Sina/Sqr Signal with Sqr  Noise * Hann * Hann
            if bitand(test_mode, 2^(curr_noise-1))
                % Mixer
                switch curr_noise
                    case {1}
                        mixer_SinaSig = sina_signal   .* sina_noise;
                        mixer_SqrSig  = square_signal .* sina_noise;
                    case {2}
                        mixer_SinaSig = sina_signal   .* sqr_noise;
                        mixer_SqrSig  = square_signal .* sqr_noise;
                    case {3}
                        mixer_SinaSig = sina_signal   .* sina_noise .* hann_win;
                        mixer_SqrSig  = square_signal .* sina_noise .* hann_win;
                    case {4}
                        mixer_SinaSig = sina_signal   .* sqr_noise .* hann_win;
                        mixer_SqrSig  = square_signal .* sqr_noise .* hann_win;
                    case {5}
                        mixer_SinaSig = sina_signal   .* sina_noise .* hann_win .* hann_win;
                        mixer_SqrSig  = square_signal .* sina_noise .* hann_win .* hann_win;
                    case {6}
                        mixer_SinaSig = sina_signal   .* sqr_noise .* hann_win .* hann_win;
                        mixer_SqrSig  = square_signal .* sqr_noise .* hann_win .* hann_win;
                end

                sina_data = sum(mixer_SinaSig);
                sqr_data  = sum(mixer_SqrSig);
                % Find max int
                if sina_data > curr_log_sina(curr_noise,1)
                    curr_log_sina(curr_noise,1) =  sina_data;
                    curr_log_sina(curr_noise,2) =  phase_step * 360 / signal_phase_cnt;
                end
                % Find min int
                if sina_data < curr_log_sina(curr_noise,3)
                    curr_log_sina(curr_noise,3) =  sina_data;
                    curr_log_sina(curr_noise,4) =  phase_step * 360 / signal_phase_cnt;
                end
                % Find max int
                if sqr_data > curr_log_sqr(curr_noise,1)
                    curr_log_sqr(curr_noise,1) =  sqr_data;
                    curr_log_sqr(curr_noise,2) =  phase_step * 360 / signal_phase_cnt;
                end
                % Find min int
                if sqr_data < curr_log_sqr(curr_noise,3)
                    curr_log_sqr(curr_noise,3) =  sqr_data;
                    curr_log_sqr(curr_noise,4) =  phase_step * 360 / signal_phase_cnt;
                end
                if int_threshold > 0 && sina_data > int_threshold
                    figure;
                    fft_signal = fft(mixer_SinaSig,NFFT);
                    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
                    plot(f/1e3, 10*log10(pxx_signal), 'b');
                    xlabel('freq[kHz]')
                    ylabel('PSD[dB]')
                    title(sprintf('SinaSig NoiseType: %d freq: %d KHz, phase: %d, INT: %d', curr_noise, noise_freq/1000, phase_step, sina_data));
                end
                if int_threshold > 0 && sqr_data > int_threshold
                    figure;
                    fft_signal = fft(mixer_SinaSig,NFFT);
                    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
                    plot(f/1e3, 10*log10(pxx_signal), 'b');
                    xlabel('freq[kHz]')
                    ylabel('PSD[dB]')
                    title(sprintf('SqrSig NoiseType: %d freq: %d KHz, phase: %d, INT: %d', curr_noise, noise_freq/1000, phase_step, sqr_data));
                end
            end
        end
    end

    test_freq(1, signal_step+1) = signal_freq;

    % Collect all measure data on current noise freq.
    for curr_noise=1:(noise_types)
        % BIT0: Sina/Sqr Signal with Sina Noise
        % BIT1: Sina/Sqr Signal with Sqr  Noise
        % BIT2: Sina/Sqr Signal with Sina Noise * Hann
        % BIT3: Sina/Sqr Signal with Sqr  Noise * Hann
        % BIT4: Sina/Sqr Signal with Sina Noise * Hann * Hann
        % BIT5: Sina/Sqr Signal with Sqr  Noise * Hann * Hann
        if bitand(test_mode, 2^(curr_noise-1))
            sina_max_log(curr_noise, signal_step+1) = curr_log_sina(curr_noise,1);
            sina_min_log(curr_noise, signal_step+1) = curr_log_sina(curr_noise,3);
            sina_maxphase_log(curr_noise, signal_step+1) = curr_log_sina(curr_noise,2);
            sina_minphase_log(curr_noise, signal_step+1) = curr_log_sina(curr_noise,4);

            sqr_max_log(curr_noise, signal_step+1) = curr_log_sqr(curr_noise,1);
            sqr_min_log(curr_noise, signal_step+1) = curr_log_sqr(curr_noise,3);
            sqr_maxphase_log(curr_noise, signal_step+1) = curr_log_sqr(curr_noise,2);
            sqr_minphase_log(curr_noise, signal_step+1) = curr_log_sqr(curr_noise,4);
        end
    end
end


%
% Fig for all [Sina Noise] Case
%
figure;
leg_array = [];
plot_cnt = 0;
for curr_noise=1:(noise_types)
    if bitand(test_mode, 2^(curr_noise-1))
        plot_cnt = plot_cnt + 1;
%        leg_array = [leg_array; noise_type_str_cell{curr_noise}];
        leg_array = [leg_array; noise_type_str(curr_noise)];
        plot(test_freq, sina_max_log(curr_noise,:), plot_color(plot_cnt));
        hold on
    end
end
if plot_cnt > 0
    legend(leg_array);
    xlabel('freq[kHz]');
    ylabel('Max Int');
    title(sprintf('Sample freq from %d Khz ~ %d Khz [step: %d] with Sina Noise', signal_freq_start/1000, signal_freq_end/1000, signal_freq_step/1000));
end


%
% Fig for all [Sqr Noise] Case
%
figure;
leg_array = [];
plot_cnt = 0;
for curr_noise=1:(noise_types)
    if bitand(test_mode, 2^(curr_noise-1))
        plot_cnt = plot_cnt + 1;
%        leg_array = [leg_array; noise_type_str_cell{curr_noise}];
        leg_array = [leg_array; noise_type_str(curr_noise)];
        plot(test_freq, sqr_max_log(curr_noise,:), plot_color(plot_cnt));
        hold on
    end
end
if plot_cnt > 0
    legend(leg_array);
    xlabel('freq[kHz]');
    ylabel('Max Int');
    title(sprintf('Sample freq from %d Khz ~ %d Khz [step: %d] with Sqr Noise', signal_freq_start/1000, signal_freq_end/1000, signal_freq_step/1000));
end


%
% Save to WS
%
assignin('base', 'signal_sweep_freq',           test_freq);
assignin('base', 'signal_sweep_sina_max',       sina_max_log);
assignin('base', 'signal_sweep_sina_max_phase', sina_maxphase_log);
assignin('base', 'signal_sweep_sqr_max',        sqr_max_log);
assignin('base', 'signal_sweep_sqr_max_phase',  sqr_maxphase_log);
assignin('base', 'signal_sweep_noise_type',     noise_type_str);

if save_ws
    disp('Save wordspace')
    file_name = sprintf('log\\%s_noise_sweep.mat', desc);
    disp(file_name);
    save(file_name, 'test_freq', 'sina_max_log', 'sina_maxphase_log', 'sqr_max_log', 'sqr_maxphase_log', 'noise_type_str');
end


%
% End
%
disp('Total sample per noise freq.');
y = n_samples;

