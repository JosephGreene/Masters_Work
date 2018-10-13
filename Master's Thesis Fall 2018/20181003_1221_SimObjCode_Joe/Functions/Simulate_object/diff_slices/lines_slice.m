% return a slice of given dimensions with two lines, one horizontal on top
% and one vertical on the right side
function slice = lines_slice(horx,hory,horlen,verx,very,verlen,...
    slice_size_px,linwidth_px)

slice = zeros(slice_size_px);

% plot horizontal lines single pixel width
for i=1:length(horx)
    slice(hory(i):hory(i)+linwidth_px,horx(i):horx(i)+horlen(i)) = 1;
end

% plot vertical lines single pixel width
for i=1:length(verx)
    slice(very(i):very(i)+verlen(i),verx(i):verx(i)+linwidth_px) = 1;
end

end