﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    interface INonLinearSolution
    {
        double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector);
    }
}
