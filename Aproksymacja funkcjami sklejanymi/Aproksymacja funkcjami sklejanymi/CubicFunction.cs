using System;

namespace Spline_Approximation_dll
{
    public class CubicFunction
    {
        public double A { get; set; }
        public double B { get; set; }
        public double C { get; set; }
        public double D { get; set; }

        public CubicFunction()
        {

        }

        public CubicFunction(double a, double b, double c, double d)
        {
            A = a;
            B = b;
            C = c;
            D = d;
        }

        public double Value(double x)
        {
            return A * x * x * x + B * x * x + C * x + D;
        }

        public double FirstDerivative(double x)
        {
            return 3 * A * x * x + 2 * B * x + C;
        }

        public static CubicFunction operator +(CubicFunction cubic, double add)
        {
            return new CubicFunction(cubic.A + add, cubic.B + add,
                                     cubic.C + add, cubic.D + add);
        }

        public static CubicFunction operator +(CubicFunction cubic, CubicFunction cubic2)
        {
            return new CubicFunction(cubic.A + cubic2.A, cubic.B + cubic2.B,
                                     cubic.C + cubic2.C, cubic.D + cubic2.D);
        }

        public static CubicFunction operator -(CubicFunction cubic, double subtract)
        {
            return new CubicFunction(cubic.A - subtract, cubic.B - subtract,
                                     cubic.C - subtract, cubic.D - subtract);
        }
        public static CubicFunction operator -(CubicFunction cubic, CubicFunction cubic2)
        {
            return new CubicFunction(cubic.A - cubic2.A, cubic.B - cubic2.B,
                                     cubic.C - cubic2.C, cubic.D - cubic2.D);
        }

        public static CubicFunction operator *(CubicFunction cubic, double multiply)
        {
            return new CubicFunction(cubic.A * multiply, cubic.B * multiply,
                                      cubic.C * multiply, cubic.D * multiply);
        }

        public static CubicFunction operator /(CubicFunction cubic, double divide)
        {
            if (divide == 0)
            {
                throw new DivideByZeroException();
            }
            return new CubicFunction(cubic.A / divide, cubic.B / divide,
                                     cubic.C / divide, cubic.D / divide);
        }

        public double[] CalculateRealRoots()
        {
            if (A != 0)
            {
                double f, g, h;
                f = (C - B * B / 3 / A) / A;
                g = ((2 * B * B / 9 / A - C) / A * B / 3 + D) / A;
                h = (g * g / 4) + (f * f * f / 27);
                if (h > 0)
                {
                    return OneRoot(g, h);

                }
                else if (f == 0 && g == 0)
                {
                    return TripleRoot();
                }
                else
                {
                    return ThreeRoots(g, h);
                }
            }
            return new double[0];

        }
        private double[] OneRoot(double g, double h)
        {
            double x1 = -(g / 2) + Math.Sqrt(h);
            double x2 = -(g / 2) - Math.Sqrt(h);
            double x3 = B / (3 * A);
            double root1 = Math.Cbrt(x1) + Math.Cbrt(x2) - x3;
            double[] roots = new double[1];
            roots[0] = root1;
            return roots;
        }

        private double[] TripleRoot()
        {
            double[] roots = { Math.Cbrt(D / A) };
            return roots;
        }

        private double[] ThreeRoots(double g, double h)
        {
            double i = Math.Sqrt((g * g / 4) - h);
            double j = Math.Cbrt(i);
            double k = Math.Acos(-g / (2 * i));
            double m = Math.Cos(k / 3);
            double n = Math.Sqrt(3) * Math.Sin(k / 3);
            double p = (-B) / (3 * A);
            double root1 = (2 * j * m) + p;
            double root2 = (-j * (m + n)) + p;
            double root3 = (-j * (m - n)) + p;
            double[] roots = new double[3];
            roots[0] = root1;
            roots[1] = root2;
            roots[2] = root3;
            return roots;
        }
    }
}
