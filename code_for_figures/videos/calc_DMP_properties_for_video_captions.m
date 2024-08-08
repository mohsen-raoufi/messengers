clc;
clear all;

p2expltArr = exp(-[1:2:100]./5);
p2msngrArr = p2expltArr(1:end-1);

i_p2e_list = [ 1, 21, 39, 49, 15, 07];
i_p2m_list = [49, 30, 49, 01, 02, 12];

ind = 6;
i_p2e = i_p2e_list(ind);
i_p2m = i_p2m_list(ind);


p2e = p2expltArr(i_p2e);
p2m = p2msngrArr(i_p2m);

disp(["p2e: ", p2e, ", p2m: ", p2m])
disp(["log10 p2e: ", log10(p2e), ", log10 p2m: ", log10(p2m)])

m = p2m/(p2e+p2m);
disp(["m: ", m])

tau_e = inf;
tau_m = inf;
div_zero = true;
if(p2e~=0)
    tau_e = 1/p2e;
    div_zero = false;
end
if(p2m~=0)
    tau_m = 1/p2m;
    div_zero = false;
end

disp(["tau_e: ", tau_e, ", tau_m: ", tau_m])
disp(["log10 tau_e: ", log10(tau_e), ", log10 tau_m: ", log10(tau_m)])

if(~div_zero)
    tau_s = 0.5*(p2m+p2e)./(p2e.*p2m);
    disp(["tau_s: ", tau_s, ", log10 tau_s: ", log10(tau_s)])
end



% disp("")
