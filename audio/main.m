function main

[c0,fs] = audioread('Comp0.wav');
[c1,fs] = audioread('Comp1.wav');
[c2,fs] = audioread('Comp2.wav');
[c3,fs] = audioread('Comp3.wav');

% normSize = size(c0)
% size(c1)
% size(c2)
% size(c3)



% subplot(2,2,1)
plot(0:1/fs:(size(c0)-1)/fs,c0);
% spectrogram(c0(:,1),hamming(4096),256,8192,fs,'yaxis')
title('Time Domain Plot of Original Sin Sweep Signal');
xlabel('Time (s)');
ylabel('Amplitude');
saveas(gcf, '../dia/Comp0plot','pdf');

% subplot(2,2,2)
% c1 = c1(1:normSize);
plot(0:1/fs:(size(c1)-1)/fs,c1);
% spectrogram(c1(:,1),hamming(4096),256,8192,fs,'yaxis')
title('Time Domain Plot of Sin Sweep Signal with Compressor 1 applied');
xlabel('Time (s)');
ylabel('Amplitude');
saveas(gcf, '../dia/Comp1plot','pdf');


% subplot(2,2,3)
% c2 = c2(1:normSize);
plot(0:1/fs:(size(c2)-1)/fs,c2);
% spectrogram(c2(:,1),hamming(4096),256,8192,fs,'yaxis')
title('Time Domain Plot of Sin Sweep Signal with Compressor 2 applied');
xlabel('Time (s)');
ylabel('Amplitude');
saveas(gcf, '../dia/Comp2plot','pdf');


% subplot(2,2,4)
% c3 = c3(1:normSize);
plot(0:1/fs:(size(c3)-1)/fs,c3);
% spectrogram(c2(:,1),hamming(4096),256,8192,fs,'yaxis')
title('Time Domain Plot of Sin Sweep Signal with Compressor 3 applied');
xlabel('Time (s)');
ylabel('Amplitude');
saveas(gcf, '../dia/Comp3plot','pdf');


 
% plot(c1)
% pause
% plot(c2)
% pause
% plot(c3)

% spectrogram(c1(:,1),hamming(4096),256,8192,fs,'yaxis')
% pause
% spectrogram(c2(:,1),hamming(4096),256,8192,fs,'yaxis')
% pause
% spectrogram(c2(:,1),hamming(4096),256,8192,fs,'yaxis')



end
