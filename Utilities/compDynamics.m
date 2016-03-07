function compDynamics (id1, name1, id2, name2)


[cmetrics, cfourier, cgraph] = nfkbmetrics (id1);
[tmetrics, tfourier, tgraph] = nfkbmetrics (id2);

 

%Max amplitude
figure ('Name', 'Max Amplitude')
x = 0:0.25:15;
chist = histogram(cmetrics.max_amplitude,x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.max_amplitude,x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Max Amplitude --IL4 Pretreatment (2d+ IMDM Swap)');
xlabel (' Amplitude (A.U.)')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 



%Longest Consecutive Duration
figure ('Name', 'Longest Consecutive Duration')
x = 0:0.5:20;
chist = histogram(cmetrics.envelope (:,4),x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.envelope (:,4),x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Longest Consecutive Duration --IL4 Pretreatment (2d+ IMDM Swap)');
xlabel ('Hours')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 



%max_integral 
figure ('Name', 'Integral')
x = 0:2:75;
chist = histogram(cmetrics.max_integral,x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.max_integral,x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Integral --IL4 Pretreatment (2d+ IMDM Swap)');
xlabel ('Activity (A.U.)')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 

%max_derivative 
figure ('Name', 'Derivative')
x = 0:0.6:25;
chist = histogram(cmetrics.max_derivative,x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.max_derivative,x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Derivative --IL4 Pretreatment (2d+ IMDM Swap)');
xlabel ('Rate of Change of Nuclear Translocation')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 

%peak frequency
%frequency that contains most the signal power
figure ('Name', 'Peak Frequency')
x = 0:0.05:1;
chist = histogram(cmetrics.peakfreq,x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.peakfreq,x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Peak Frequency --IL4 Pretreatment (2d + IMDM Swap)');
xlabel ('Frequency (1/hr) ')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 

%Fraction oscillating

figure ('Name', 'Fraction Oscillating')
x = 0:0.05:1;
chist = histogram(cmetrics.oscfrac(:,3),x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.oscfrac(:,3),x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Fraction Oscillating--IL4 Pretreatment (2d+ IMDM Swap)');
xlabel ('Oscillatory Content')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 



%Duration above threshold
figure ('Name', 'Duration above Threshold')
x = 0:0.5:20;
chist = histogram(cmetrics.duration (:,4),x);
chist.FaceColor = 'auto';
chist.EdgeColor = 'k';
chist.Normalization = 'probability';
hold on;
thist = histogram(tmetrics.duration (:,4),x);
thist.Normalization ='probability';
thist.FaceColor = 'auto';
thist.EdgeColor='r';
title ('Duration above Threshold --IL4 Pretreatment (2d+ IMDM Swap)');
xlabel ('Hours')
ylabel ('Probability')
legend ('Ctrl (100nM CpG)', '10 ng/ml IL4 pretreatment (100nM CpG)','Location', 'northeast') 






