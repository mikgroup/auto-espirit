function noise = addNoise(sigma, shape)

  noise = normrnd(0, sigma/sqrt(2), shape) + i * normrnd(0, sigma/sqrt(2), shape);

end
