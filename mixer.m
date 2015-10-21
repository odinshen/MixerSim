%
% Mixer waveform generator
%

close all
% clc

disp('(1) Noise Sweeping, fixed Sample');
disp('(2) Sample Sweeping, fixed Sina/Sqr noise');
disp('(3) Sample Sweeping, mat noise pattern');

sim_select = input('Choise type:      ');

disp('Save MAT files ?');
disp('(1) Yes, please');
disp('(0) No, thanks');
saving = input(':');

if saving > 0
    mat_prefix = input('Input mat prefix: ');
end

disp(mat_prefix)
% return

if sim_select == 1
    disp('Type 1')
    default_setting = input('Setting (1)Default (2) Input:   ');
    if default_setting == 1
        if saving
            mixer_noise_sweeping(mat_prefix);
        else
            mixer_noise_sweeping();
        end
    elseif default_setting == 2
        test_select = input('Test mode   : ');
        mixer_freq  = input('Mixer Freq  : ');
        noise_start = input('Noise Start : ');
        noise_end   = input('Noise End   : ');
        noise_step  = input('Noise Step  : ');
        phase_cnt   = input('Phase Cnt   : ');
        noise_level = input('Noise Level : ');
        int_thres   = input('Int Thresh  : ');
        mixer_noise_sweeping(mat_prefix, test_select, mixer_freq, noise_start, noise_end, noise_step, phase_cnt, noise_level, int_thres);
    else
        if saving
            mixer_noise_sweeping(mat_prefix, 63,200e3,1e3,1e6,200,36,1,0)  %% quick test
        else
            mixer_noise_sweeping('quick', 63,200e3,1e3,1e6,200,36,1,0)  %% quick test
        end
    end
elseif sim_select == 2
    disp('Type 2')
    default_setting = input('Setting (1)Default (2) Input:   ');
    if default_setting == 1
        if saving
            mixer_signal_sweeping(mat_prefix);
        else
            mixer_signal_sweeping()
        end
    elseif default_setting == 2
        test_select = input('Test mode    : ');
        noise_freq  = input('Noise Freq   : ');
        signal_start = input('Signal Start : ');
        signal_end   = input('Signal End   : ');
        signal_step  = input('Signal Step  : ');
        phase_cnt   = input('Phase Cnt    : ');
        noise_level = input('Noise Level  : ');
        int_thres   = input('Int Thresh   : ');
        mixer_signal_sweeping(mat_prefix, test_select, noise_freq, signal_start, signal_end, signal_step, phase_cnt, noise_level, int_thres);
    else
        if saving
            mixer_signal_sweeping(mat_prefix, 63, 32.75e3, 1e3, 1e6, 200, 36, 1, 0)  %% quick test
        else
            mixer_signal_sweeping('quick', 63, 32.75e3, 1e3, 1e6, 200, 36, 1, 0)  %% quick test
        end
    end
elseif sim_select == 3
    disp('Type 3')
    default_setting = input('Setting (1)Default (2) Input:   ');
    if default_setting == 1
        if saving
            mixer_signal_sweeping_500kmat(mat_prefix);
        else
            mixer_signal_sweeping_500kmat();
        end
    elseif default_setting == 2
        test_select  = input('Test mode    : ');
        noise_scale  = input('Noise Scale  : ');
        signal_start = input('Signal Start : ');
        signal_end   = input('Signal End   : ');
        signal_step  = input('Signal Step  : ');
        phase_cnt    = input('Phase Cnt    : ');
        noise_level  = input('Noise Level  : ');
        int_thres    = input('Int Thresh   : ');
        mixer_signal_sweeping_500kmat(mat_prefix, test_select, noise_scale, signal_start, signal_end, signal_step, phase_cnt, noise_level, int_thres);
    else
        if saving
            mixer_signal_sweeping_500kmat(mat_prefix, 3*16+8, 1, 10e3, 20e3, 200, 4, 1, 0)  %% quick test
        else
            mixer_signal_sweeping_500kmat('quick', 3*16+8, 1, 10e3, 11e3, 200, 4, 1, 0)  %% quick test
        end
    end

elseif sim_select == 9
    if saving
        disp('Type Test 1')
        mixer_noise_sweeping(mat_prefix, 63,200e3,1e3,1e6,200,36,1,0)  %% quick test
        disp('Type Test 2')
        mixer_signal_sweeping(mat_prefix, 63, 32.75e3, 1e3, 1e6, 200, 36, 1, 0)  %% quick test
        disp('Type Test 3')
        mixer_signal_sweeping_500kmat(mat_prefix, 3*16+8, 1, 10e3, 100e3, 1e3, 4, 1, 0)  %% quick test
    else
        disp('Type Test 1')
        mixer_noise_sweeping('quick', 63,200e3,1e3,1e6,200,36,1,0)  %% quick test
        disp('Type Test 2')
        mixer_signal_sweeping('quick', 63, 32.75e3, 1e3, 1e6, 200, 36, 1, 0)  %% quick test
        disp('Type Test 3')
        mixer_signal_sweeping_500kmat('quick', 3*16+8, 1, 10e3, 100e3, 1e3, 4, 1, 0)  %% quick test
    end
end

