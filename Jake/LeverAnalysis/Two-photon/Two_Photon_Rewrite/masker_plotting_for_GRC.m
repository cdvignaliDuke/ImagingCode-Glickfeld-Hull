% make a mask for the neurons in the example 2P figure for the last figure
% of the GRC 2017 poster

poster_mask = mask3D(:,:,[39, 53, 60]);
figure; subplot(3,1,1); suptitle('img94');
imagesc(squeeze(poster_mask(:,:,1))); truesize;
title('neuron 39');
subplot(3,1,2);
imagesc(squeeze(poster_mask(:,:,2))); truesize;
 title('neuron 53');
subplot(3,1,3);
imagesc(squeeze(poster_mask(:,:,3))); truesize;
 title('neuron 60');

poster_mask = sum(poster_mask, 3);
figure; imagesc(poster_mask);
truesize
title('img94 example neurons 39, 53, 60')