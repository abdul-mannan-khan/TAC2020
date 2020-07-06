function vout = enlarge(vin,pct)
    c = mean(vin);
    Delta = diff(vin);
    vout = c+[-Delta Delta]*pct/2;
end