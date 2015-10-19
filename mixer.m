%
% Mixer waveform generator
%

close all
% clc

disp('(1)Noise Sweeping, fixed Sample');
disp('(2)Sample Sweeping, fixed Sina/Sqr noise');
disp('(3)Sample Sweeping, mat noise pattern');

sim_select = input('Choise type:    ');

if sim_select == 1
    disp('Type 1')
    default_setting = input('Setting (1)Default (2) Input:   ');
    if default_setting == 1
        mixer_noise_sweeping()
    elseif default_setting == 2
        test_select = input('Test mode   : ');
        mixer_freq  = input('Mixer Freq  : ');
        noise_start = input('Noise Start : ');
        noise_end   = input('Noise End   : ');
        noise_step  = input('Noise Step  : ');
        phase_cnt   = input('Phase Cnt   : ');
        noise_level = input('Noise Level : ');
        int_thres   = input('Int Thresh  : ');
        mixer_noise_sweeping(test_select, mixer_freq, noise_start, noise_end, noise_step, phase_cnt, noise_level, int_thres);
        %{
            test_mode        = (2^noise_types) - 1;    % Input 1, Select all bits
            mixer_freq       = 200e3; % 200K           % Input 2
            noise_freq_start = 1e3;   % 1K             % Input 3
            noise_freq_end   = 2e6;   % 2M            % Input 4
            noise_freq_step  = 200;   % 200 % 1e3 1k   % Input 5
            noise_phase_cnt  = 36;    % 360/18         % Input 6
            noise_level      = 1;
            int_threshold    = 0;
        %}
    else
        mixer_noise_sweeping(3,200e3,1e3,40e3,200,36,1,0)  %% quick test
    end
elseif sim_select == 2
    disp('Type 2')
    default_setting = input('Setting (1)Default (2) Input:   ');
    if default_setting == 1
        mixer_signal_sweeping()
    elseif default_setting == 2
        test_select = input('Test mode    : ');
        noise_freq  = input('Noise Freq   : ');
        signal_start = input('Signal Start : ');
        signal_end   = input('Signal End   : ');
        signal_step  = input('Signal Step  : ');
        phase_cnt   = input('Phase Cnt    : ');
        noise_level = input('Noise Level  : ');
        int_thres   = input('Int Thresh   : ');
        mixer_signal_sweeping(test_select, noise_freq, signal_start, signal_end, signal_step, phase_cnt, noise_level, int_thres);
        %{
            test_mode         = (2^noise_types) - 1;    % Input 1, Select all bits
            noise_freq        = 32.75e3; % 32.75K       % Input 2
            signal_freq_start = 1e3;     % 1K           % Input 3
            signal_freq_end   = 2e6;     % 2M           % Input 4
            signal_freq_step  = 200;   % 200 % 1e3 1k   % Input 5
            signal_phase_cnt  = 36;    % 360/18         % Input 6
            noise_level       = 1;                      % Input 7
            int_threshold     = 0;                      % Input 8
        %}
    else
        mixer_signal_sweeping(63, 32.75e3, 10e3, 500e3, 200, 36, 1, 0)  %% quick test
    end
elseif sim_select == 3
    disp('Type 3')
    default_setting = input('Setting (1)Default (2) Input:   ');
    if default_setting == 1
        mixer_signal_sweeping_500kmat()
    elseif default_setting == 2
        test_select  = input('Test mode    : ');
        noise_scale  = input('Noise Scale  : ');
        signal_start = input('Signal Start : ');
        signal_end   = input('Signal End   : ');
        signal_step  = input('Signal Step  : ');
        phase_cnt    = input('Phase Cnt    : ');
        noise_level  = input('Noise Level  : ');
        int_thres    = input('Int Thresh   : ');
        mixer_signal_sweeping_500kmat(test_select, noise_freq, signal_start, signal_end, signal_step, phase_cnt, noise_level, int_thres);
        %{
            test_mode         = (2^noise_types) - 1;    % Input 1, Select all bits
            noise_freq        = 32.75e3; % 32.75K       % Input 2
            signal_freq_start = 1e3;     % 1K           % Input 3
            signal_freq_end   = 2e6;     % 2M           % Input 4
            signal_freq_step  = 200;   % 200 % 1e3 1k   % Input 5
            signal_phase_cnt  = 36;    % 360/18         % Input 6
            noise_level       = 1;                      % Input 7
            int_threshold     = 0;                      % Input 8
        %}
    else
        mixer_signal_sweeping_500kmat(16+8, 1, 200e3, 1e6, 500, 36, 1, 0)  %% quick test
    end

elseif sim_select == 9
    disp('Type Test 1')
    mixer_noise_sweeping(3,200e3,1e3,40e3,200,36,1,0)  %% quick test
    disp('Type Test 2')
    mixer_signal_sweeping(3, 32.75e3, 1e3, 40e3, 200, 36, 1, 0)  %% quick test
    disp('Type Test 3')
    mixer_signal_sweeping_500kmat(3, 1, 200e3, 1e6, 500, 6, 1, 0)  %% quick test
end

