from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr

def solve_P2():
  session = WolframLanguageSession()
  session.evaluate("""
Clear["Global`*"];
M = 4866; m = 2433; k2 = 80000; c2 = 10000; omega = 2.2143; k1 = 1025 * 9.8 * Pi;
f = 4890; c1 = 167.8395; m0 = 1165.992;
start = 510; stop = 520; dim = 1;
P2[cm2_] := Sum[cm2 * ((x1'[i] - x2'[i])^2) * dim, {i, start, stop, dim}] /. 
  NDSolve[{m*x2''[t]==k2*(x1[t]-x2[t])+cm2*(x1'[t]-x2'[t]), 
    f*Cos[omega * t]==k2*(x1[t]-x2[t])+cm2*(x1'[t]-x2'[t])+k1*x1[t]+c1*x1'[t]+m0*x1''[t]+M*x1''[t], 
    x1[0]==x2[0]==x1'[0]==x2'[0]==0},
  {x1, x2}, {t, start, stop}]
""")
  
  
  def PP2(x, d):
    return session.evaluate(wl.Evaluate(wlexpr(f'First[P2[{x}]-P2[{x-d}]]')))
  
  def find(start, stop, f, esp=2*1e-5, n=100):
    d = (stop - start) / n
    print(f"find({start}, {stop}, {d})")
    m = (start + stop) / 2
    v = f(m, d)
    if abs(v) < esp:
      return [m, v]
    if v < 0:
      return find(start, m, f)
    return find(m, stop, f)

  return find(0, 100000, PP2)

# solve_P2()
"""
find(0, 100000, 1000.0)
find(0, 50000.0, 500.0)
find(25000.0, 50000.0, 250.0)
find(25000.0, 37500.0, 125.0)
find(31250.0, 37500.0, 62.5)
find(34375.0, 37500.0, 31.25)
find(34375.0, 35937.5, 15.625)
find(34375.0, 35156.25, 7.8125)
...
find(34473.31529058459, 34473.315290585306, 7.130438461899757e-12)
find(34473.31529058459, 34473.31529058495, 3.5652192309498786e-12)
"""
# [34473.31529058477, 0.0]