bifdxl = load('attr_pool.txt');
%bifdxr = load('bd-x-rl.txt');

hold on
plot(bifdxl(:,1),bifdxl(:,2),'.','Color',[0/255,0/255,255/255],'MarkerSize',4,'DisplayName','bifd')
%plot(bifdxr(:,1),bifdxr(:,2),'.','Color',[0/255,191/255,255/255],'MarkerSize',4,'DisplayName','bifdxr')

legend('off')
