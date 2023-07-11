A = readtable("SCOPE_01.CSV");
B = readtable("SCOPE_02.CSV");
C = readtable("SCOPE_03.CSV");

%Now let us do some processing
figure 
plot(A{1:1400, 2} - min(A{1:1400, 2}))
hold
plot((A{1:1400, 3} - min(A{1:1400, 3})) / (max(A{1:1400, 3}) - min(A{1:1400, 3})))
figure 
plot(B{1:2000, 2} - min(B{1:2000, 2}))
hold
plot((B{1:2000, 3} - min(B{1:2000, 3})) / (max(B{1:2000, 3}) - min(B{1:2000, 3})))
figure 
plot(C{1:900, 2} - min(C{1:900, 2}))
hold
plot((C{1:900, 3} - min(C{1:900, 3})) / (max(C{1:900, 3}) - min(C{1:900, 3})))

%Auto Correlation
x1 = A{:, 3};
x2 = B{:, 3};
x3 = C{:, 3};
x1 = (x1 - min(x1));%/(max(x1) - min(x1));
x2 = (x2 - min(x2));%/(max(x2) - min(x2));
x3 = (x3 - min(x3));%/(max(x3) - min(x3));

[y1, lags1] = xcorr(x1);% autocorr(x1, NumLags=3999);
[y2, lags2] = xcorr(x2);% autocorr(x2, NumLags=3999);
[y3, lags3] = xcorr(x3);% autocorr(x3, NumLags=3999);

figure
plot(lags1(:), y1)
figure
plot(lags2, y2)
figure
plot(lags3, y3)

%x^6
pnseq1 = comm.PNSequence('Polynomial',[6 1 0], 'InitialConditions', [1 1 1 1 1 1], 'SamplesPerFrame', 63 * 4);
pnseq2 = comm.PNSequence('Polynomial',[6 5 4 3 2 1 0], 'InitialConditions', [1 1 1 1 1 1], 'SamplesPerFrame', 63 * 4);

out1 = pnseq1();
out2 = pnseq2();

[y4, lags4] = xcorr(out1);
[y5, lags5] = xcorr(out2);

figure
plot(lags4, y4)
figure
plot(lags5, y5)