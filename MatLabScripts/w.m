w1 = load('w1.txt');
w2 = load('w2.txt');
w3 = load('w3.txt');

hold on
plot(w1(:,1),w1(:,2),'LineWidth',2,'DisplayName','w1');
plot(w2(:,1),w2(:,2),'LineWidth',2,'DisplayName','w2');
plot(w3(:,1),w3(:,2),'LineWidth',2,'DisplayName','w3');

legend('off')