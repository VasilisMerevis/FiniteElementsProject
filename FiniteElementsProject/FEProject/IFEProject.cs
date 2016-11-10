using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    interface IFEProject
    {
        /// <summary>
        /// Changes the default name of .txt input file
        /// </summary>
        /// <param name="newInputFile"></param>
        void ChangeInputFile(string newInputFile);
        void ReadInputFile();
        void CreateModel();
    }
}
