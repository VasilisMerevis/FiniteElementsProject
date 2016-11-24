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

            
        }
    }
}
