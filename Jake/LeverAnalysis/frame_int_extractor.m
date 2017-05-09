function extracted_frames = frame_int_extractor(full_movie, index_of_frames, num_frames_before_indexed, num_frames_after_indexed);
%looks at the entire movie and extracts only those frames which lie within
%a specified window surrounding each indexed frame. So if frame #10 is
%indexed then it will collect frames 5:12 if num_frames_before_indexed=5
%and num_frames_after_indexed=2
    extracted_frames = [];
    for frame_num = 1:length(index_of_frames); 
        extracted_frames = cat(3, extracted_frames, full_movie(:,:, [index_of_frames(frame_num)+num_frames_before_indexed:index_of_frames(frame_num)+num_frames_after_indexed]));
    end
end