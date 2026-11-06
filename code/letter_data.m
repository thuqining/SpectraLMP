function [x, dt, meta] = letter_data(csvPath)
% LETTER_DATA  Case-study net load (v23): CoSSMic industrial 1 and
% residential 5 and 6, with measured rooftop PV at industrial 1 and
% residential 6, no OPSD gap-filled samples. The first 60 days train the
% policies and the following 365 days test them (425 days in total, from
% 2016-07-28, inside the 441-day common gap-free window).
[x, dt, meta] = load_data(csvPath, 425/30, {'ind1', 'res5', 'res6'});
end
