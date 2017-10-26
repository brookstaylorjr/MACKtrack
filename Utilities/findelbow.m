function elbow = findelbow(x,y)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% elbow = findelbow(x,y)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FINDELBOW tries to fit this shape to a vector: _/ - the best-fit "elbow point" will be returned
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



v1 = 3:length(x)-3;
tot_err = zeros(size(v1));
for j = 1:length(v1)
    x1 = x(1:v1(j)); y1 = y(1:v1(j));
    f1 = polyfit(x1,y1,1);
    err1 = sum(abs(polyval(f1,x1)-y1));
    x2 = x(v1(j):end); y2 = y(v1(j):end);
    f2 = polyfit(x2,y2,1);
    err2 = sum(abs(polyval(f2,x2)-y2));
    tot_err(j) = err1+err2;
end
elbow = x(v1(find(tot_err==min(tot_err),1,'last')));
