
dest_dir = 'c:\GoogleDrive\syncUpFolder\odinMatlabCode\mixer_Sim\'

%{
file_list = ['mixer.m'; 'mixer_noise_sweeping.m'; 'mixer_signal_sweeping.m'; 'mixer_signal_sweeping_500kmat.m'; 'sync_to_gdrive.m'];
noise_type = ['mixer.m'; 'mixer_noise_sweeping.m'; 'mixer_signal_sweeping.m'; 'mixer_signal_sweeping_500kmat.m'; 'sync_to_gdrive.m';
              'LCD_noise_0kohm.mat'; 'LCD_noise_20kohm.mat'; 'LCD_noise_High_Z.mat';
              'square_noise_0kohm_fs50MHz.mat'; 'square_noise_20kohm_fs50MHz.mat'; 'square_noise_High_Z_fs50MHz.mat'];
files = size(file_list);
%}

status = zeros(1,20);
[status(1,1),mess,messid] = copyfile('mixer.m', dest_dir);
[status(1,2),mess,messid] = copyfile('mixer_noise_sweeping.m', dest_dir);
[status(1,3),mess,messid] = copyfile('mixer_signal_sweeping.m', dest_dir);
[status(1,4),mess,messid] = copyfile('mixer_signal_sweeping_500kmat.m', dest_dir);
[status(1,5),mess,messid] = copyfile('sync_to_gdrive.m', dest_dir);

% [status(1,6),mess,messid] = copyfile('LCD_noise_0kohm.mat', dest_dir);
% [status(1,7),mess,messid] = copyfile('LCD_noise_20kohm.mat', dest_dir);
% [status(1,8),mess,messid] = copyfile('LCD_noise_High_Z.mat', dest_dir);
% [status(1,9),mess,messid] = copyfile('square_noise_0kohm_fs50MHz.mat', dest_dir);
% [status(1,10),mess,messid] = copyfile('square_noise_20kohm_fs50MHz.mat', dest_dir);
% [status(1,11),mess,messid] = copyfile('square_noise_High_Z_fs50MHz.mat', dest_dir);

if sum(status) == 5
    disp('Copy Succ~~');
else
    check_file = find(status<1, 1);
    disp('Something Wrong on file: ');
    disp(check_file);
end


