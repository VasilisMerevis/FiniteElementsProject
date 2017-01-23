using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public class FEProject : IFEProject
    {
        private InputData data;

        public FEProject()
        {
            this.data = new InputData();
        }

        /// <summary>
        /// Changes the default name of .txt input file
        /// </summary>
        /// <param name="newInputFile"></param>
        public void ChangeInputFile(string newInputFile)
        {
            data.SetNewInputFile(newInputFile);
        }

        public void ReadInputFile()
        {
            data.ReadAllData();
        }

        public void CreateModel()
        {
            if (data.elementType.Contains("NLBeam") | data.elementType.Contains("NLTruss"))
            {
                Console.WriteLine("Non-linear Solution is:");
                IAssembly Exercise1Frame = new Assembly(data);
                Exercise1Frame.GetStiffnessMatrices();
                Exercise1Frame.BoundedDOFsVector = data.boundaryDof;

                ISolver nonLin = new StaticSolver();
                nonLin.AssemblyData = Exercise1Frame;
                nonLin.LinearScheme = new CholeskyFactorization();
                //nonLin.SetSolutionMethodToCholesky();
                nonLin.SetNonLinearMethodToLoadControlledNewtonRaphson();
                
                nonLin.ActivateNonLinearSolver = true;
                nonLin.SolveStatic(data.externalForcesVector);
                nonLin.PrintSolution();
            }
           
            else
            {
                IAssembly Exercise1Frame = new Assembly(data);
                Exercise1Frame.BoundedDOFsVector = data.boundaryDof;
                Exercise1Frame.GetStiffnessMatrices();
                Exercise1Frame.InitializeMatrices();

                Console.WriteLine("Linear Solution is:");
                ISolver newSolu = new StaticSolver();
                //newSolu.SetSolutionMethodToPCG();
                newSolu.LinearScheme = new PCGSolver();
             
                newSolu.AssemblyData = Exercise1Frame;
                newSolu.SolveStatic(data.externalForcesVector);
                newSolu.PrintSolution();
            }
        }
      
    }
}
