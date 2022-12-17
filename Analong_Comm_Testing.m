x = 0:0.1:10;
y = x.^2;
figure; plot(x, y); grid on;
saveas(gcf, 'test_graph.png');
