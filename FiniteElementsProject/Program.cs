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
            FEProject linearCantilever = new FEProject();
            linearCantilever.ReadInputFile();
            linearCantilever.CreateModel();

            FEProject nLCantilever = new FEProject();
            nLCantilever.ChangeInputFile("Writelines2NLCantilever.txt.");
            nLCantilever.ReadInputFile();
            nLCantilever.CreateModel();
        }
    }
}
