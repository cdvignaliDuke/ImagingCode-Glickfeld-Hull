function K = Toeplitz_matrix (ker, Nframes);
K=zeros(Nframes);
for i=1:Nframes
    range=i:min(Nframes,i+length(ker)-1);
    K(i,range) = ker(1:length(range));
end