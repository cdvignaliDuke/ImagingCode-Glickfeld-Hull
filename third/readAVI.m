
aviObj = VideoReader('001_000_000_eye.avi');

nt = aviObj.Duration * aviObj.FrameRate; % get #frame

mov = zeros(aviObj.Height, aviObj.Width, nt, 'uint8'); % initialize to uint8

for k = 1:nt
    mov(:,:,k) = rgb2gray(readFrame(aviObj)); % convert rgb to gray
end