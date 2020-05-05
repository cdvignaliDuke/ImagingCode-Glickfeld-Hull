
function [allclusters, allclusters_init, centroidcorr, dendmem, dunnsinitial]=meta_k_means(eventsmat, distmetric,sc,nclust)

if nargin==1;
    distmetric='correlation';
    sc = 0;
    nclust = 5;
end

expr = eventsmat;
numden = size(expr,1);
numtimes = size(expr, 2);
counts = zeros(numden); 
centroidmat = [];


% k clusters
if sc==0
    n =nclust;
else
    n=size(sc,1);
end
for k = n;%4
    
    % repeat t times
    for t = 1:1000
        if rem(t,100) == 0
            fprintf([num2str(t) '\n'])
        end
        
%         if t == 1
%         pool = parpool;                      % Invokes workers
%         stream = RandStream('mlfg6331_64');  % Random number stream
%         options = statset('UseParallel',1,'UseSubstreams',1,...
%             'Streams',stream);
%         end
        if sc == 0
            [IDX, ~, ~, ~] = kmeans(expr, k,  'Distance', distmetric);
        else
            [IDX, ~, ~, ~] = kmeans(expr, k, 'Distance', distmetric, 'Start', sc);
        end
%         [IDX, ~, ~, ~] = kmeans(expr, k, 'Options', options, 'Distance', distmetric);
        
        % update counts
        % returns a matrix showing the number of times each pair of
        % dendrites were clustered together
       
        for i = 1:numden
            counts(i,IDX == IDX(i)) = counts(i,IDX == IDX(i)) + 1;
        end
        
    end
    
    if ~isempty(IDX);
        
        thr = 0.8; % 800 out of 1000 times
        
%         cooccur_matrix = counts/t;
%         
%         clusters = [];
%         for xx = 1:numden
%             for yy = 1:xx-1
%                 if counts(xx, yy) >= thr*t
%                     clusters = [clusters xx];
%                     clusters = [clusters yy];
%                 end
%             end
%         end
%         
%         % sort list of dendrites
%         clusters = sort(clusters);
%         
%         % delete the repeated dendrites from the list
%         list = unique(clusters);
        [xx, yy] = find(counts >= thr*t);
        leftC = yy;
        cnt = 1;
        while ~isempty(leftC)
            kk = xx(yy == leftC(1));
            leftC = setdiff(leftC,kk);
            allclusters{cnt, 1} = kk';
            allclusters{cnt,2} = sum(expr(kk,:), 1)/length(kk); % centroid
            
            cnt = cnt + 1;
        end
        
        
        allclusters_init = allclusters;
        % vector of average point to centroid distances for each cluster
%         numvect = [];
        
 
        %             centroidcorr = corrcoef(centroidmat');
        
        combclusters=0;
        
        [~, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, combclusters);
        
        endclustering=0;
        
        flag = 0;
        rep = 1;
        
        while endclustering==0;
            combineclusters=[];
            numclusters=size(allclusters,1);
            fprintf(['Combine rep = ' num2str(rep) '- dunnsinitial = ' num2str(dunnsinitial) '- ' num2str(numclusters) 'clusters \n'])
            for i=1:numclusters-1
                for j=(i+1):numclusters;
                    if j > size(allclusters,1)
%                         combineclusters = [];
                        break
                    end
                    [dunns]=ioclustervalidity2_min(allclusters, expr, [i j]);
                    
                    if dunns >= dunnsinitial;
                        combineclusters=[i j];
                        dunnsinitial=dunns;
                        fprintf(['[' num2str(i) ' ' num2str(j) '] - dunnsinitial = ' num2str(dunnsinitial) '\n'])
                    end
                end
                
            end
            rep = rep+1;
            if ~isempty(combineclusters);
                [allclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, [combineclusters(1) combineclusters(2)]);
            else
                endclustering=1;
            end
            if numclusters <= k || size(allclusters,1) <= k;
                endclustering=1;
            end
        end
        
        %             numclusters=size(allclusters,1);
        
        %             for zz=1:numclusters;
        %                 clusters=allclusters{zz,1};
        %                 for zzz=1:length(clusters);
        %                     clusters(zzz)=clusters(zzz)+sum(silentdendcnt(1:clusters(zzz)));
        %                 end
        %                 allclusters{zz,1}=clusters;
        %             end
        
    else
        allclusters=[];
        centroidcorr=[];
        dendmem=[];
        dunnsinitial=0;
    end
    
end



