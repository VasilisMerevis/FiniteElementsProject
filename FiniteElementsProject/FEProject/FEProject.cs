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
            if (data.elementType.Contains("NLBeam"))// | data.elementType.Contains("NLTruss"))
            {
                //Creation of local, global and total stiffness matrices
                Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                Exercise1Frame.GetStiffnessMatrices();
                //Exercise1Frame.CreateTotalStiffnessMatrix();
                NLSolver solution = new NLSolver(Exercise1Frame, 1000, data);
                solution.SolveWithMethod("Load Controlled Newton-Raphson");
                VectorOperations.PrintVector(solution.solutionVector);
            }
            else if (data.elementType.Contains("NLTruss"))
            {
                //Creation of local, global and total stiffness matrices
                Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                Exercise1Frame.GetStiffnessMatrices();

                //Exercise1Frame.CreateTotalStiffnessMatrix();
                NLSolver solution = new NLSolver(Exercise1Frame, 10000, data);
                solution.SolveWithMethod("Newton-Raphson");
                VectorOperations.PrintVector(solution.solutionVector);
            }
            else
            {
                //Creation of local, global and total stiffness matrices
                Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                Exercise1Frame.GetStiffnessMatrices();
                Exercise1Frame.InitializeMatrices();
                Exercise1Frame.CreateTotalStiffnessMatrix();
                Exercise1Frame.GetMassMatrices();
                //Creation reduced matrix depended on boundary conditions
                double[,] reducedTotalStiff = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame.TotalStiffnessMatrix, data.boundaryDof);

                //Solution using Cholesky factorization with forward and backward substitution
                DirectSolver solution = new DirectSolver(reducedTotalStiff, data.externalForcesVector);
                solution.SolveWithMethod("Cholesky");
                Console.WriteLine();

                VectorOperations.PrintVector(solution.GetSolutionVector);

                Console.WriteLine();
                IterativeSolver solution2 = new IterativeSolver(reducedTotalStiff, data.externalForcesVector, 1000);
                solution2.SolveWithMethod("PCG");

                VectorOperations.PrintVector(solution2.GetSolutionVector);




                Console.WriteLine();
                Discretization2DFrame Exercise1Frame2 = new Discretization2DFrame(data);
                Exercise1Frame2.GetStiffnessMatrices();
                Exercise1Frame2.InitializeMatrices();
                Exercise1Frame2.CreateTotalStiffnessMatrix();
                Exercise1Frame2.GetMassMatrices();
                //Creation reduced matrix depended on boundary conditions
                double[,] reducedTotalStiff2 = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame2.TotalStiffnessMatrix, data.boundaryDof);

                LinearSolver solution3 = new DirectMethods();
                solution3.SetSolutionMethodToGauss();
                //solution3.GetStiffnessMatrixAndForceVector(reducedTotalStiff, data.externalForcesVector);
                
                solution3.Solve(reducedTotalStiff2, data.externalForcesVector);
                solution3.PrintSolution();

                //VectorOperations.PrintVector(solution3.GetSolutionVector);
                Console.WriteLine("newSolu is:");
                StaticSolver newSolu = new StaticSolver();
                
                //ILinearSolution newSolu = new LinearSolution();
                newSolu.SetSolutionMethodToGauss();
                newSolu.Solve(reducedTotalStiff2, data.externalForcesVector);
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
