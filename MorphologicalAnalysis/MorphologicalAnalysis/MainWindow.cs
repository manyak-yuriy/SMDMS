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

            initP[0] = V.DenseOfArray(new [] {0.3, 0.5, 0.2});
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

            int permInd = 0;
            int currFIndex = N;


        }
    }
}
