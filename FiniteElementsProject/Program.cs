using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
    class Program
    {
        static void Main(string[] args)
        {
            IFEProject linearCantilever = new FEProject();
            linearCantilever.ReadInputFile();
            linearCantilever.CreateModel();

            IFEProject nLCantilever = new FEProject();
            nLCantilever.ChangeInputFile("Writelines2NLCantilever.txt.");
            nLCantilever.ReadInputFile();
            nLCantilever.CreateModel();

            IFEProject linearTruss = new FEProject();
            linearTruss.ChangeInputFile("LinearTruss.txt");
            linearTruss.ReadInputFile();
            linearTruss.CreateModel();

            double[] initialU = { 0, 0 };
            double[] initialDotU = { 0, 0 };
            double[] initialDdotU = { 0, 10 };
            double initialTime = 0;
            double totalTime = 2.8;
            int timeStepsNumber = 10;
            double[,] Kmatrix = { { 6, -2 }, { -2, 4 } };
            double[,] Mmatrix = { { 2, 0 }, { 0, 1 } };
            double[] Fvector = { 0, 10 };
            CentralDifferencesSolver dynamic = new CentralDifferencesSolver(initialTime, initialU, initialDotU, initialDdotU, totalTime, timeStepsNumber, Kmatrix, Mmatrix, Fvector);
            dynamic.SolveExplicit();
            Console.WriteLine();
            dynamic.PrintExplicitSolution();

          
        }
    }
}
