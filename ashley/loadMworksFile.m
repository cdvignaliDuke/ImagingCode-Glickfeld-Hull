function mWorksData = loadMworksFile(subjectNumber,exptDate,exptTime,varargin)
if ~isempty(varargin)
    mWorksDataPath = varargin{1};
else
    mWorksDataPath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\Behavior\Data';

end
mworksFileName = ['data-' 'i' subjectNumber '-' exptDate '-' exptTime]; 
mWorksDataPlusBackup = load(fullfile(mWorksDataPath,mworksFileName));
mWorksData = mWorksDataPlusBackup.input;
end