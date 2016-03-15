function newton(f, g!, H!, x0)
  x = copy(x0); n = length(x)
  gx = zeros(n); Hx = zeros(n,n)
  g!(x, gx) # gx = g(x)
  H!(x, Hx) # Hx = H(x)
  while norm(gx) > 1e-8
    BLAS.axpy!(-1, Hx\gx, x)
    g!(x, gx) # gx = g(x)
    H!(x, Hx) # Hx = H(x)
  end
  return x
end
