function maxIndex=syncFMCWSymbol(y,offsetIndex,y0,K)

L=1*K;
corr=zeros(1,L);
maxIndex = 0;
for i=1+offsetIndex:L+offsetIndex
    if(length(y) < K-1 + i)
        maxIndex = -1
        break;
    end
    corr(i-offsetIndex)=y0*y((0:K-1)+i)';
end
if (maxIndex  == 0)
    figure;  plot(abs(corr));
    abs_corr = abs(corr);
    [m, idx] = max(abs_corr);
    [sorted_abs_corr, sorted_idx] = sort(abs_corr, 'descend');
    maxIndex = idx + offsetIndex;
end
%maxIndex = sorted_idx(2) + offsetIndex; % try second max
% maxCorr=0;
% maxIndex=0;
% for i=1:L
%     if abs(corr(i))>maxCorr
%         maxCorr=abs(corr(i));
%         maxIndex=i+offsetIndex;
%     end
% end