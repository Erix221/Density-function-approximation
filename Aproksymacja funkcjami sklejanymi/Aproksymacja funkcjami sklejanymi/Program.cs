using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Numerics;


namespace Aproksymacja_funkcjami_sklejanymi
{
    class Program
    {
        static void Main(string[] args)
        {           
            SplineApproximation sp = new SplineApproximation();
            double[,] Testresult = sp.TestSplineSinus();

            Console.WriteLine("Porównanie wartości");
            for (int i = 0; i < 100; i++)
            {
                Console.WriteLine("Wartość aproksymowana: " + Testresult[i, 0] + "| Wartość dokładna: " + Testresult[i, 1]);
                Console.WriteLine("Różnica: " + (Testresult[i, 1] - Testresult[i, 0]));
            }

        }
    }
}
