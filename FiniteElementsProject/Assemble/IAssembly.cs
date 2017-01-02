using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public interface IAssembly
    {
        void InitializeMatrices();
        Element1D[] GetStiffnessMatrices();
        void GetMassMatrices();
        void UpdateValues(double[] totalDisplacementVector);
        double[,] CreateTotalStiffnessMatrix();
        double[,] CreateTotalMassMatrix();
        double[] CreateTotalInternalForcesVector();
        int[] BoundedDOFsVector
        {
            get;
            set;
        }
    }
}
