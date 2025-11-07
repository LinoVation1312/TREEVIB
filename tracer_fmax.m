load("C:\Users\lconord.LYD\OneDrive - MCO\Documents\Peuplier\fmax_results.mat")

smoothing = 50;

fmax1est_raw = export_data.Arbre1.Arbre1_Est.fmax;
fmax1nord_raw = export_data.Arbre1.Arbre1_Nord.fmax;
fmax2est_raw = export_data.Arbre2.Arbre2_Est.fmax;
fmax2nord_raw = export_data.Arbre2.Arbre2_Nord.fmax;

fmax1est = smoothdata(fmax1est_raw, 'movmean', smoothing);
fmax1nord = smoothdata(fmax1nord_raw, 'movmean', smoothing);
fmax2est = smoothdata(fmax2est_raw, 'movmean', smoothing);
fmax2nord = smoothdata(fmax2nord_raw, 'movmean', smoothing);


dates = export_data.Arbre1.Arbre1_Est.dates; 

if isnumeric(dates)
    dates = datetime(dates, 'ConvertFrom', 'datenum');
elseif ischar(dates) || isstring(dates)
    dates = datetime(dates, 'InputFormat', 'dd-MMM-yyyy');
end

figure;
plot(dates, fmax1est, 'b-', 'LineWidth', 1.5);
hold on;
plot(dates, fmax1nord, 'r-', 'LineWidth', 1.5);
plot(dates, fmax2est, 'g-', 'LineWidth', 1.5);
plot(dates, fmax2nord, 'm-', 'LineWidth', 1.5);

xlabel('Date');
ylabel('fmax (lissé)');
title(['fmax en fonction de la date']);
legend('Arbre1 Est', 'Arbre1 Nord', 'Arbre2 Est', 'Arbre2 Nord');
grid on;

datetick('x', 'dd-mmm-yyyy', 'keepticks');