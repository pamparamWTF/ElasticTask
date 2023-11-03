using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static ElasticTask.Elements;

namespace ElasticTask
{
    class SLAE
    {
        public int N;
        public int[] ig;
        public int[] jg;
        public double[] di;
        public double[] ggl;
        public double[] ggu;
        public double[] F;

        // for Solver
        private double[] l, u, d;
        private double[] temp, temp0;
        private double[] r, z, p;

        public double[] q { get; set; }
        public SLAE(Elements elements, Mesh mesh)
        {
            N = 3 * mesh.BigMesh[0].isMesh.Count * mesh.BigMesh[1].isMesh.Count * mesh.BigMesh[2].isMesh.Count;
            Portrait(elements);

            int size = ig[N];

            l = new double[size];
            u = new double[size];
            d = new double[N];
            q = new double[N];
            temp = new double[N];
            temp0 = new double[N];
            r = new double[N];
            z = new double[N];
            p = new double[N];

            elements.BuildGlobalMatrix(ig, jg, ref di, ref ggl);

            SetBoundarySecond(elements);
            SetBoundaryFirst(elements);
        }
        private void Portrait(Elements elements)
        {
            List<SortedSet<int>> map = new List<SortedSet<int>>();

            for (int i = 0; i < N; i++)
            {
                map.Add(new SortedSet<int>());
            }

            foreach (Element elem in elements.elements)
                for (int i = 0; i < elem.Node_global.Count; i++)
                    for (int j = 0; j < elem.Node_global.Count; j++)
                        if (i > j)
                            map[elem.Node_global[i]].Add(elem.Node_global[j]);

            ig = new int[map.Count + 1];

            ig[0] = 0;

            for (int i = 0; i < map.Count; i++)
            {
                ig[i + 1] = ig[i] + map[i].Count;
            }

            jg = new int[ig[ig.Length - 1]];

            for (int i = 0; i < map.Count; i++)
            {
                var jind = map[i].ToArray();
                for (int j = 0; j < jind.Length; j++)
                    jg[ig[i] + j] = jind[j];
            }

            di = new double[N];
            ggl = new double[jg.Length];
            ggu = new double[jg.Length];
            F = new double[N];

            for (int i = 0; i < N; i++)
            {
                di[i] = 0;
                F[i] = 0;
            }

            for (int i = 0; i < jg.Length; i++)
            {
                ggl[i] = 0;
                ggu[i] = 0;
            }

        }
        private void SetBoundaryFirst(Elements elements)
        {

        }
        private void SetBoundarySecond(Elements elements)
        {
            double integral, hx, hy, hz;
            foreach (BoundaryConditions bc in elements.boundaryConditions)
            {
                if (bc.type == 2)
                {
                    foreach (Element element in bc.elements2D)
                    {
                        hx = element.hx();
                        hy = element.hy();
                        hz = element.hz();
                        for (int i = 0; i < 4; i++)
                        { 
                            integral = element.BoundaryIntegral(i, bc.Ngran, hx, hy, hz);

                            foreach (int node in element.Node_global)
                            {
                                //F[node] += 
                            }
                        }
                    }
                }
            }

        }

    }
}
