function newton(f, g, H, x)
  while norm(g(x)) > 1e-8
    x = x - H(x)\g(x)
  end
  return x
end
