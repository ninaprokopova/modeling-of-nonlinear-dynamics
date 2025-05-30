%data = load('..\data\phase_p.txt'); DisplayName = 'phase_p'; type = ".-";
data = load('..\data\ellipse.txt'); DisplayName = 'ellipse'; type = ".";
%data = load('..\data\lp.txt'); DisplayName = 'lp'; type = ".-";
%data = load('..\data\rot_num.txt'); DisplayName = 'rot_num'; type = ".-";
%data = load('..\data\attrp_mu_x.txt'); DisplayName = 'attrp_mu_x'; type = ".";
%data = load('..\data\bifd.txt'); DisplayName = 'bifd'; type = ".";

hold on
plot(data(:,1), data(:,2), ...
    type, ...
    'MarkerSize', 6, ...
    'Color', [77/255 190/255 238/255], ...
    'DisplayName', DisplayName)

legend('off')
