double DoubleFactorial(int a) {
  double b=1;
  while(a>1) {
    b=b*a;
    a-=2;
  }
  return b;
}
