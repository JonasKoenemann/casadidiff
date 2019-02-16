
int main()
{
  a = adlib.sym();
  b = adlib.sym();
  f = adlib.add(adlib.mul(a,b),a)

  f.eval_fcn({a,b},{4,2});

  F = Function({a,b},{f});
  F(4,2);

}
