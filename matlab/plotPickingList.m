points = load('../data/picking_list.txt');
plot3(points(:,1), points(:,2), points(:,3));
hold on;
pts = load('../data/Heli_position.txt');
plot3(pts(:,1), pts(:,2), pts(:,3));
axis equal;