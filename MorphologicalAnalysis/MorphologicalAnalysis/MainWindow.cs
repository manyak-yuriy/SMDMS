using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace MorphologicalAnalysis
{
    public partial class MainWindow : Form
    {


        public MainWindow()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            var M = Matrix<double>.Build;
            var V = Vector<double>.Build;

            int N = 3;

            // n[i] is a number of alternatives for i-th factor
            int[] n = new int[N];

            n[0] = 3;
            n[1] = 4;
            n[2] = 2;

            // i-th element is a vector of initial probabilities for i-the factor's alternatives
            Vector<double>[] initP = new Vector<double>[N];

            initP[0] = V.DenseOfArray(new[] { 0.3, 0.5, 0.2 });
            initP[1] = V.DenseOfArray(new[] { 0.4, 0.3, 0.1, 0.2 });
            initP[2] = V.DenseOfArray(new[] { 0.3, 0.7 });

            int n_max = n.Max();

            var cMatrix = new Matrix<double>[N, N];

            for (int i = 0; i < N - 1; i++)
                for (int j = i + 1; j < N; j++)
                    cMatrix[i, j] = M.Sparse(n[i], n[j]);

            // F1-F2 block
            cMatrix[0, 1][0, 0] = 0.5;
            cMatrix[0, 1][2, 1] = -0.5;

            // F1-F3 block
            cMatrix[0, 2][0, 0] = 0.2;
            cMatrix[0, 2][1, 0] = 0.3;

            // F2-F3 block
            cMatrix[1, 2][0, 0] = 0.5;
            cMatrix[1, 2][2, 1] = -1;




            // F[i] permutations


            int totalPermCnt = 1;

            foreach (int num in n)
                totalPermCnt *= num;

            // create F[i] to store permutations

            int[][] F = new int[N][];

            for (int i = 0; i < N; i++)
            {
                F[i] = new int[totalPermCnt];
                F[i][0] = 0;
            }


            int permInd = 1;
            int currFIndex = N - 1;

            do
            {

                for (int i = 0; i < N; i++)
                    F[i][permInd] = F[i][permInd - 1];

                while (++F[currFIndex][permInd] == n[currFIndex])
                {
                    F[currFIndex][permInd] = 0;
                    currFIndex--;
                }

                currFIndex = N - 1;

                ++permInd;
            }
            while (permInd < totalPermCnt);

            // -- all permutations generated 

            permInd = 0;

            // calculate column of c values
            double[] c = new double[totalPermCnt];

            while (permInd < totalPermCnt)
            {
                c[permInd] = 1;

                for (int m = 0; m < N - 1; m++)
                    for (int l = m + 1; l < N; l++)
                        c[permInd] *= (cMatrix[m, l][F[m][permInd], F[l][permInd]] + 1);

                permInd++;
            }
            // --

            double[,] pMarginal = new double[N, totalPermCnt];

            for (int i = 0; i < N; i++)
                for (permInd = 0; permInd < totalPermCnt; permInd++)
                {
                    pMarginal[i, permInd] = c[permInd];

                    for (int Find = 0; Find < N; Find++)
                        if (Find != i)
                            pMarginal[i, permInd] *= initP[Find][F[Find][permInd]];
                }
            // -- calculated marginal probs

            // calc normalized marginal probs inplace
            // normlize only inside current alternatives' values

            for (int i = 0; i < N; i++)
            {
                for (int currFValue = 0; currFValue < n[i]; currFValue++)
                {
                    double sum = 0;
                    for (permInd = 0; permInd < totalPermCnt; permInd++)
                        if (currFValue == F[i][permInd])
                            sum += pMarginal[i, permInd];

                    for (permInd = 0; permInd < totalPermCnt; permInd++)
                        if (currFValue == F[i][permInd])
                            pMarginal[i, permInd] /= sum;
                }
            }

            // calculate P[i] matrices

            Matrix<double>[] PMatrices = new Matrix<double>[N];

            for (int currInd = 0; currInd < N; currInd++)
            {

                int nextInd = currInd + 1;

                if (nextInd == N)
                    nextInd = 0;

                PMatrices[currInd] = Matrix<double>.Build.Dense(n[currInd], n[nextInd]);

                for (int k = 0; k < n[currInd]; k++)
                    for (int l = 0; l < n[nextInd]; l++)
                    {
                        PMatrices[currInd][k, l] = 0;

                        permInd = 0;

                        while (permInd < totalPermCnt)
                        {
                            if (F[currInd][permInd] == k && F[nextInd][permInd] == l)
                                PMatrices[currInd][k, l] += pMarginal[nextInd, permInd];
                            permInd++;
                        }
                    }

                // -- filled the PMatrices[currInd] 
            }
            // -- filled all the PMatrices


            // calculate the matrices' product

            Matrix<double> P = Matrix<double>.Build.DenseIdentity(n[0], n[0]);

            for (int i = 0; i < N; i++)
                P = P * PMatrices[i];

            Vector<double> x = Vector<double>.Build.Dense(n[0]);
            x[0] = 1;

            Vector<double> xNext = x;


            // find x1 iteratively
            do
            {
                x = xNext;
                xNext = P * x;
            } while ((x - xNext).L2Norm() > 1e-6);

            Vector<double>[] xArray = new Vector<double>[N];

            xArray[N - 1] = PMatrices[N - 1] * x;

            for (int i = N - 2; i >= 0; i--)
                xArray[i] = PMatrices[i] * xArray[i + 1];

            // Normalize - apparently, unnecessary

            //for (int i = 0; i < N; i++)
            //    xArray[i] /= xArray[i].L1Norm();

            for (int i = 0; i < N; i++)

        }
    }
}
