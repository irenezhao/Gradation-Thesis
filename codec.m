function cod = codec(m)

g1 = [1 1 1];
g2 = [1 0 1];

m1 = conv(m,g1);
m1 = mod(m1,2);

m2 = conv(m,g2);

m2 = mod(m2,2);

l = length(m1);
cod = zeros(1,2*l);

for i=1:l;
    cod(2*i-1) = m1(i);
    cod(2*i) = m2(i);
end