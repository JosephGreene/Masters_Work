% this fuction returns a slice of given size with the text inserted in the
% slice at the given position
function slice = text_slice(text,slice_size_px,xpos,ypos)
text_im = text2im(text);
[n,m] = size(text_im);
slice = zeros(slice_size_px);
for i = 1:length(xpos)
    slice(ypos(i):ypos(i)+(n-1),xpos(i):xpos(i)+(m-1)) = not(text_im);
end
end