function result = getMeasures(h,test_label)

s = size(h,1);
result = [];
for i=1:s
    f = h(i,:);
    result = [result; performanceMeasure(test_label, f)];
end



