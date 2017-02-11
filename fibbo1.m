function v=fibbo1(n)
F= [1 1];
raksha (n)


function raksha(n)
if numel(F)<n
    raksha(n-1)
    F(n)=F(n-1)+F(n-2);
end
end
v=F(n);
end
   