% this fuction returns a slice of given size with the text inserted in the
% slice at the given position
function slice = text_sliceasd(text,slice_size_px,xpos,ypos)
text_im = text2im(text);
[m,n] = size(text_im);
slice = zeros(slice_size_px);
slice(xpos:xpos+(m-1),ypos:ypos+(n-1)) = not(text_im);
end