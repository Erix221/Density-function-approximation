using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace Aproksymacja_funkcjami_sklejanymi
{
    public class Tests
    {

        public double[][] irisdataset()
        {
            double[][] iris = new double[2][];
            iris[0] = new double[46];
            iris[0][0] = 5.1;
            iris[0][1] = 4.9;
            iris[0][2] = 4.7;
            iris[0][3] = 4.6;
            iris[0][4] = 5;
            iris[0][5] = 5.4;
            iris[0][6] = 4.6;
            iris[0][7] = 5;
            iris[0][8] = 4.4;
            iris[0][9] = 4.9;
            iris[0][10] = 5.4;
            iris[0][11] = 4.8;
            iris[0][12] = 4.8;
            iris[0][13] = 4.3;
            iris[0][14] = 5.8;
            iris[0][15] = 5.7;
            iris[0][16] = 5.4;
            iris[0][17] = 5.1;
            iris[0][18] = 5.4;
            iris[0][19] = 5.1;
            iris[0][20] = 4.6;
            iris[0][21] = 5.1;
            iris[0][22] = 4.8;
            iris[0][23] = 5;
            iris[0][24] = 5;
            iris[0][25] = 5.2;
            iris[0][26] = 5.2;
            iris[0][27] = 4.7;
            iris[0][28] = 4.8;
            iris[0][29] = 5.4;
            iris[0][30] = 5.2;
            iris[0][31] = 5.5;
            iris[0][32] = 4.9;
            iris[0][33] = 5;
            iris[0][34] = 4.4;
            iris[0][35] = 5.1;
            iris[0][36] = 5;
            iris[0][37] = 4.5;
            iris[0][38] = 4.4;
            iris[0][39] = 5;
            iris[0][40] = 5.1;
            iris[0][41] = 4.8;
            iris[0][42] = 5.1;
            iris[0][43] = 4.6;
            iris[0][44] = 5.3;
            iris[0][45] = 5;

            return iris;
        }

        public double[,] TestSplineSinus(double[] x)
        {
            double[] y = new double[1000000];
            for (int i = 0; i < 1000000; i++) // wyliczam wartości Sinusa
            {
                y[i] = Math.Sin(x[i]);
            }
            SplineApproximation sp = new SplineApproximation(0, 2 * Math.PI, 6, x, y); // zmieniac trzecią wartość, jest to numer przedzialów
            double h = (sp.EndOfInterval - sp.StartOfInterval) / sp.NumberOfIntervals;
            var sw1 = Stopwatch.StartNew();
            double res1 = 0;
            double[] result1 = sp.Approximate1();
            Console.WriteLine("New: " + sw1.ElapsedMilliseconds + " ms:  " + res1);
            double[,] sinsprox = new double[1000000, 2];
            for (int i = 0; i < 1000000; i++)
            {
                sinsprox[i, 0] = sp.FunctionValue1(result1, x[i]); // szybsze rozwiazanie
                sinsprox[i, 1] = y[i];
            }
            return sinsprox;
        }

        public double[,] TestSplineSinus1(double[] x)
        {
            double[] y = new double[1000000];
            for (int i = 0; i < 1000000; i++) // wyliczam wartości Sinusa
            {
                y[i] = Math.Sin(x[i]);
            }
            SplineApproximation sp = new SplineApproximation(0, 2 * Math.PI, 6, x, y); // zmieniac trzecią wartość, jest to numer przedzialów
            double h = (sp.EndOfInterval - sp.StartOfInterval) / sp.NumberOfIntervals;
            var sw1 = Stopwatch.StartNew();
            double res1 = 0;
            double[] result1 = sp.Approximate1();
            Console.WriteLine("New: " + sw1.ElapsedMilliseconds + " ms:  " + res1);
            double[,] sinsprox = new double[1000000, 2];
            for (int i = 0; i < 1000000; i++)
            {
                sinsprox[i, 0] = sp.FunctionValue1(result1, x[i]); // szybsze rozwiazanie
                sinsprox[i, 1] = y[i];
            }
            return sinsprox;
        }

        public double[][] DensityFunctionTest()
        {
            Random r = new Random();
            double range = 5;
            double rDouble = r.NextDouble() * range; //for doubles
            double[][] x = irisdataset();
            DensitiyFunction df = new DensitiyFunction(4, 6, 0.2, 0.1, x[0]);
            double[][] result = df.CalculateBinary();
            return result;
        }

        public double[][] DensityFunctionTest1()
        {
            Random r = new Random();
            double range = 5;
            double rDouble = r.NextDouble() * range; //for doubles
            double[][] x = irisdataset();
            DensitiyFunction df = new DensitiyFunction(4, 6, 0.2, 0.1, x[0]);
            double[][] result = df.CalculateSmart();
            return result;
        }
    }
    }
