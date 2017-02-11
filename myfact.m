function v = myfact(n)
    s = ['myfact('n')'];
    disp(s);
    a=1;
    if(n>0)
        s = ['enteing if for n='n];
        disp (s);
        a = n * myfact(n-1);
        s = ['exiting if for n ='n];
        disp(s);
    end
    v =a;
end

