function y = getDescriptors(linea, block, mb_size, theta, set2, set3, H, xEnter, yEnter, xExit, yExit)

diff = variance(linea, block, mb_size, theta);

y = [theta; H; diff; mean(set2(:)); mean(set3(:)); var(set2(:)); var(set3(:)); xEnter; yEnter; xExit; yExit];

end