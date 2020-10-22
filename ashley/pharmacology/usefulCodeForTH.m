%% use during 2P imaging
quickRet.m %quick analysis of 2P retinotopy experiments
quickImgFromSbx.m %save a tiff of the first 100 frames of a scanbox run

%% data organization
DART_V1_PV_contrast

%% data processing and multiday
regdata_selectcells
multiday_findsamecells

%% response analysis
spontaneousActivity_multiday
contrastResp_nocellmatch_multiday
oriTuning_nocellmatch_multiday
