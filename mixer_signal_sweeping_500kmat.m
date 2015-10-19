%
% Mixer waveform generator
% Sweep [SIGNAL] from noise_freq_start ~ noise_freq_end
% Noise from mat file
%

function y = mixer(test_mode, noise_scale, signal_freq_start, signal_freq_end, signal_freq_step, signal_phase_cnt, noise_level, int_threshold)
%function y = mixer(test_mode, mixer_freq, noise_freq_start, noise_freq_end, noise_freq_step, noise_phase_cnt, noise_level, int_threshold)

% clear
close all
% clc

fs = 50e6 * 50;   % 50M
% mixing_duration = 200e-6;   % 200us
% mixing_duration = 10e-3;   % 1ms
mixing_duration = 200e-6;   % 200us

noise_type = ['LCD 0kohm '; 'LCD 20kohm'; 'LCD High Z'; 'SQR 0kohm '; 'SQR 20kohm'; 'SQR High Z'];
plot_color = ['b'; 'r'; 'c'; 'k'; 'g'; 'm'];

noise_types = size(noise_type);
noise_types = noise_types(1);
            % BIT0: Sina/Sqr Signal with LCD_noise_0kohm
            % BIT1: Sina/Sqr Signal with LCD_noise_20kohm
            % BIT2: Sina/Sqr Signal with LCD_noise_High_Z
            % BIT3: Sina/Sqr Signal with square_noise_0kohm_fs50MHz
            % BIT4: Sina/Sqr Signal with square_noise_20kohm_fs50MHz
            % BIT5: Sina/Sqr Signal with square_noise_High_Z_fs50MHz
measure_types = 4;
            % MAX_INT
            % MAX_INT_PHASE
            % MIN_INT
            % MIN_INT_PHASE

noise_type_str = cellstr(noise_type);

if nargin < 8
    disp('[Usage]:')
    disp('        Input1: test_mode')
    disp('        Input2: noise_scale')
    disp('        Input3: sinal_freq_start')
    disp('        Input4: sinal_freq_end')
    disp('        Input5: sinal_freq_step')
    disp('        Input6: sinal_phase_cnt')
    disp('        Input7: noise_level')
    disp('        Input8: threshold')
    % Default setting
    test_mode         = (2^noise_types) - 1;   % Input 1, Select all bits
%    test_mode         = 1+2+4;                  % Input 1, Select BIT0 BIT1 BIT2
    noise_scale       = 1;        % scale by mat file      % Input 2
    signal_freq_start = 50e3;     % 1K          % Input 3
%    signal_freq_end   = 2e6;    % 2M          % Input 4
    signal_freq_end   = 400e3;   % 400K        % Input 4
    signal_freq_step  = 1e3;   % 200 % 1e3 1k   % Input 5
    signal_phase_cnt  = 6;    % 360/18         % Input 6
    noise_level       = 1;
    int_threshold     = 0;
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
disp(sprintf('Noise scale %d', noise_scale));
disp(sprintf('Sweep Sample from %d Khz ~ %d Khz [step: %d Hz] [phase_cnt: %d]', signal_freq_start/1000, signal_freq_end/1000, signal_freq_step, signal_phase_cnt));
disp(sprintf('      noise level %d', noise_level));
disp(sprintf('      check INT Threshold %d', int_threshold));

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
hann_win      = hann(n_samples)';
f             = linspace(0,fs/2,NFFT/2);
% plot(t, hann_win, 'b');
hann_hann_win = hann_win .* hann_win;


%
% Load mat waveform
%
t2  = linspace(1,n_samples,n_samples);

%mat_noise1 = [];
%mat_noise2 = [];
%mat_noise3 = [];
%mat_noise4 = [];
%mat_noise5 = [];
%mat_noise6 = [];

if bitand(test_mode, (2^0 + 2^1 + 2^2))
    leg_array = [];
    plot_cnt = 0;
    figure;
    if bitand(test_mode, 2^(0))
        plot_cnt = plot_cnt + 1;
        load('LCD_noise_0kohm.mat');
        mat_noise1 = LCD_noise_0kohm(1:n_samples)';
        plot(t2, mat_noise1, plot_color(plot_cnt));
        leg_array = [leg_array; 'LCD 0kohm '];
        hold on
    end
    if bitand(test_mode, 2^(1))
        plot_cnt = plot_cnt + 1;
        load('LCD_noise_20kohm.mat');
        mat_noise2 = LCD_noise_20kohm(1:n_samples)';
        plot(t2, mat_noise2, plot_color(plot_cnt));
        leg_array = [leg_array; 'LCD 20kohm'];
        hold on
    end
    if bitand(test_mode, 2^(2))
        plot_cnt = plot_cnt + 1;
        load('LCD_noise_High_Z.mat');
        mat_noise3 = LCD_noise_High_Z(1:n_samples)';
        plot(t2, mat_noise3, plot_color(plot_cnt));
        leg_array = [leg_array; 'LCD High Z'];
        hold on
    end
    xlabel('Time');
    ylabel('Amp');
    legend(leg_array);
    title('LCD Noise waveform');
end
if bitand(test_mode, 2^(0))
    figure;
    fft_signal = fft(mat_noise1,NFFT);
    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
    plot(f/1e3, 10*log10(pxx_signal), 'b');
    xlabel('freq[kHz]')
    ylabel('PSD[dB]')
    title('Mat Noise1 FFT');
end
if bitand(test_mode, 2^(1))
    figure;
    fft_signal = fft(mat_noise2,NFFT);
    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
    plot(f/1e3, 10*log10(pxx_signal), 'b');
    xlabel('freq[kHz]')
    ylabel('PSD[dB]')
    title('Mat Noise2 FFT');
end
if bitand(test_mode, 2^(2))
    figure;
    fft_signal = fft(mat_noise3,NFFT);
    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
    plot(f/1e3, 10*log10(pxx_signal), 'b');
    xlabel('freq[kHz]')
    ylabel('PSD[dB]')
    title('Mat Noise3 FFT');
end

if bitand(test_mode, (2^3 + 2^4 + 2^5))
    leg_array = [];
    plot_cnt = 0;
    figure;
    if bitand(test_mode, 2^(3))
        plot_cnt = plot_cnt + 1;
        load('square_noise_0kohm_fs50MHz.mat');
        mat_noise4 = square_noise_0kohm_fs50MHz(1:n_samples)';
        plot(t2, mat_noise4, plot_color(plot_cnt));
        leg_array = [leg_array; 'SQR 0kohm '];
        hold on
    end
    if bitand(test_mode, 2^(4))
        plot_cnt = plot_cnt + 1;
        load('square_noise_20kohm_fs50MHz.mat');
        mat_noise5 = square_noise_20kohm_fs50MHz(1:n_samples)';
        plot(t2, mat_noise5, plot_color(plot_cnt));
        leg_array = [leg_array; 'SQR 20kohm'];
        hold on
    end
    if bitand(test_mode, 2^(5))
        plot_cnt = plot_cnt + 1;
        load('square_noise_High_Z_fs50MHz.mat');
        mat_noise6 = square_noise_High_Z_fs50MHz(1:n_samples)';
        plot(t2, mat_noise6, plot_color(plot_cnt));
        leg_array = [leg_array; 'SQR High Z'];
        hold on
    end
    xlabel('Time');
    ylabel('Amp');
    legend(leg_array);
    title('Sqr Noise waveform');
end

if bitand(test_mode, 2^(3))
    figure;
    fft_signal = fft(mat_noise4,NFFT);
    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
    plot(f/1e3, 10*log10(pxx_signal), 'b');
    xlabel('freq[kHz]')
    ylabel('PSD[dB]')
    title('Mat Noise4 FFT');
end
if bitand(test_mode, 2^(4))
    figure;
    fft_signal = fft(mat_noise5,NFFT);
    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
    plot(f/1e3, 10*log10(pxx_signal), 'b');
    xlabel('freq[kHz]')
    ylabel('PSD[dB]')
    title('Mat Noise5 FFT');
end
if bitand(test_mode, 2^(5))
    figure;
    fft_signal = fft(mat_noise6,NFFT);
    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
    plot(f/1e3, 10*log10(pxx_signal), 'b');
    xlabel('freq[kHz]')
    ylabel('PSD[dB]')
    title('Mat Noise6 FFT');
end


for signal_step=0:signal_freq_cnt
% for noise_step=0:2
    signal_freq = signal_freq_start + signal_freq_step * signal_step;

    if rem(signal_freq, 10e3) == 0
        disp(sprintf('Sweep %d Khz', signal_freq/1000)); %10KHz
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
                        mixer_SinaSig = sina_signal   .* mat_noise1;
                        mixer_SqrSig  = square_signal .* mat_noise1;
                    case {2}
                        mixer_SinaSig = sina_signal   .* mat_noise2;
                        mixer_SqrSig  = square_signal .* mat_noise2;
                    case {3}
                        mixer_SinaSig = sina_signal   .* mat_noise3;
                        mixer_SqrSig  = square_signal .* mat_noise3;
                    case {4}
                        mixer_SinaSig = sina_signal   .* mat_noise4;
                        mixer_SqrSig  = square_signal .* mat_noise4;
                    case {5}
                        mixer_SinaSig = sina_signal   .* mat_noise5;
                        mixer_SqrSig  = square_signal .* mat_noise5;
                    case {6}
                        mixer_SinaSig = sina_signal   .* mat_noise6;
                        mixer_SqrSig  = square_signal .* mat_noise6;
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
                    title(sprintf('SinaSig NoiseType: %d noise scale: %d, phase: %d, INT: %d', curr_noise, noise_freq_scale, phase_step, sina_data));
                end
                if int_threshold > 0 && sqr_data > int_threshold
                    figure;
                    fft_signal = fft(mixer_SinaSig,NFFT);
                    pxx_signal = abs(fft_signal(1:NFFT/2)).^2;
                    plot(f/1e3, 10*log10(pxx_signal), 'b');
                    xlabel('freq[kHz]')
                    ylabel('PSD[dB]')
                    title(sprintf('SqrSig NoiseType: %d noise scale: %d, phase: %d, INT: %d', curr_noise, noise_freq_scale, phase_step, sqr_data));
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
% Fig for all [Sina Signal] Case
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
    title(sprintf('Sina Sample from %d Khz ~ %d Khz [step: %d] on different Noise Pattern', signal_freq_start/1000, signal_freq_end/1000, signal_freq_step/1000));
end


%
% Fig for all [Sqr Signal] Case
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
    title(sprintf('Sqr Sample from %d Khz ~ %d Khz [step: %d] on different Noise pattern', signal_freq_start/1000, signal_freq_end/1000, signal_freq_step/1000));
end


%
% Save to WS
%
assignin('base', 'mat_signal_sweep_noise_freq',     test_freq);
assignin('base', 'mat_signal_sweep_sina_max',       sina_max_log);
assignin('base', 'mat_signal_sweep_sina_max_phase', sina_maxphase_log);
assignin('base', 'mat_signal_sweep_sqr_max',        sqr_max_log);
assignin('base', 'mat_signal_sweep_sqr_max_phase',  sqr_maxphase_log);
assignin('base', 'mat_signal_sweep_noise_type',     noise_type_str);


%
% End
%
disp('Total sample per noise freq.');
y = n_samples;

