function [  ] = MakeVideoAVI( frames, frameRate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

video = VideoWriter('Video_1dim_GPUCB.avi' );
video.FrameRate = frameRate;
open(video);
writeVideo(video, frames);
close(video);

end

