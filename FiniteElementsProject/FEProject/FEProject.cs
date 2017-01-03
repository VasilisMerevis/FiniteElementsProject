using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public class FEProject : IFEProject
    {
        //Declaration of basic frame info
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
                nonLin.SetSolutionMethodToCholesky();
                nonLin.SetNonLinearMethodToLoadControlledNewtonRaphson(Exercise1Frame);
                nonLin.NLSolve(data.externalForcesVector);
                nonLin.PrintSolution();
            }
           
            else
            {
                ////Creation of local, global and total stiffness matrices
                //Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                //Exercise1Frame.GetStiffnessMatrices();
                //Exercise1Frame.InitializeMatrices();
                //Exercise1Frame.CreateTotalStiffnessMatrix();
                //Exercise1Frame.GetMassMatrices();
                ////Creation reduced matrix depended on boundary conditions
                //double[,] reducedTotalStiff = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame.TotalStiffnessMatrix, data.boundaryDof);

                IAssembly Exercise1Frame = new Assembly(data);
                Exercise1Frame.BoundedDOFsVector = data.boundaryDof;
                Exercise1Frame.GetStiffnessMatrices();
                Exercise1Frame.InitializeMatrices();

                double[,] reducedTotalStiff = Exercise1Frame.CreateTotalStiffnessMatrix();
                Exercise1Frame.GetMassMatrices();

                Console.WriteLine("Linear Solution is:");
                ISolver newSolu = new StaticSolver();
                newSolu.SetSolutionMethodToPCG();
                newSolu.Solve(reducedTotalStiff, data.externalForcesVector);
                newSolu.PrintSolution();
            }
        }

        


            

            

        //    //Creation of local, global and total stiffness matrices
        //    Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
        //Exercise1Frame.GetStiffnessMatrices();
        //    Exercise1Frame.CreateTotalStiffnessMatrix();
        //    Exercise1Frame.GetMassMatrices();
        //    Exercise1Frame.CreateTotalMassMatrix();
        //    //Creation reduced matrix depended on boundary conditions
        //    double[,] reducedTotalStiff = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame.TotalStiffnessMatrix, data.boundaryDof);
        //double[,] reducedMassMatrix = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame.TotalMassMatrix, data.boundaryDof);

        //double[] reducedInitialU = new double[12];
        //double[] reducedInitialV = new double[12]; reducedInitialV[9] = -1.3888;
        //    double[] reducedInitialA = new double[12];
        //CentralDifferencesSolver expSolu = new CentralDifferencesSolver(0, reducedInitialU, reducedInitialV, reducedInitialA, 1, 10, reducedTotalStiff, reducedMassMatrix);

        //expSolu.SolveExplicit();
        //    VectorOperations.PrintVector(expSolu.GetExplicitSolution);
        //    //Solution using Cholesky factorization with forward and backward substitution




            //Console.WriteLine();
            //IterativeSolver solution2 = new IterativeSolver(reducedTotalStiff, data.externalForcesVector, 1000);
            //solution2.SolveWithMethod("PCG");

            //VectorOperations.PrintVector(solution2.GetSolutionVector);
            //Console.WriteLine();

            //Matrix<double> bill = new Matrix<double>(new double[,] { { 1, 3 }, { 4, 6 } });
            //bill.PrintMatrix();
            //double a = 2;
            //Matrix<double> bill2 = a * bill;
            //Console.WriteLine();
            //bill2.PrintMatrix();
    }
}
