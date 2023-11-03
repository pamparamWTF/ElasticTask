using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Xml.Linq;

namespace ElasticTask
{
    public class BoundaryCoefs
    {
        public int Ngran { get; set; }
        public int Type { get; set; }
        public double Px { get; set; }
        public double Py { get; set; }
        public double Pz { get; set; }

        public BoundaryCoefs(int Ngran, int Type, double Px, double Py, double Pz)
        {
            this.Ngran = Ngran;
            this.Type = Type;
            this.Px = Px;
            this.Py = Py;
            this.Pz = Pz;
        }
    }

    public class GeometryData
    {
        public GeometryData(int N) 
        {
            this.N = N;
            NameGeometry = "Geometry" + N.ToString();
            SizeX = 1; SizeY = 1; SizeZ = 1;
            NameMaterial = "NewMaterial" + N.ToString();
            Ex = 1;
            Ey = 1;
            Ez = 1;
            Vxy = 1;
            Vxz = 1;
            Vyz = 1;
            Gxy = 1;
            Gxz = 1;
            Gyz = 1;
            OnModel = false;
            isotropie = true;
        }
        public int N { get; set; }
        public string NameGeometry { get; set; }
        public double SizeX { get; set; }
        public double SizeY { get; set; }
        public double SizeZ { get; set; }
        public string NameMaterial { get; set; }
        //Young's module
        public double Ex { get; set; }
        public double Ey { get; set; }
        public double Ez { get; set; }
        //Poisson's ratios
        public double Vxy { get; set; }
        public double Vxz { get; set; }
        public double Vyz { get; set; }
        //Shear modulus
        public double Gxy { get; set; }
        public double Gxz { get; set; }
        public double Gyz { get; set; }
        public bool OnModel { get; set; }
        public bool isotropie { get; set; }
        public double[,] Dmatrix { get; set; }
        public void BuildD()
        {
            Dmatrix = CalcD();
        }
        public double[,] CalcD()
        {
            double[,] matrix = new double [6, 6];

            if (isotropie)
            {
                double E = Ex, 
                    V = Vxy, 
                    G = E / (2 * (1 + V)), 
                    h = 1 - 3 * V * V - 2 * V * V * V,
                    diagD = (E / h) * (1 - V * V),
                    ggluD = (E / h) * (V + V * V);

                for (int i = 0; i < 6; i++)
                {
                    if (i < 3) matrix[i, i] = diagD;
                    else matrix[i, i] = G;
                    for (int j = 0; j < 6; j++)
                        if (i < 3 && j < 3 && i!=j) matrix[i, j] = ggluD;
                        else if (i != j) matrix[i, j] = 0;
                }
            }
            return matrix;
        }
    }
    public class MeshCoefs
    {
        public MeshCoefs(double hX, double hY, double hZ, double kX, double kY, double kZ)
        {
            this.hX = hX;
            this.hY = hY;
            this.hZ = hZ;
            this.kX = kX;
            this.kY = kY;
            this.kZ = kZ;
        }
        public double hX { get; set; }
        public double hY { get; set; }
        public double hZ { get; set; }
        public double kX { get; set; }
        public double kY { get; set; }
        public double kZ { get; set; }
    }
    public class Elements
    {
        public List<Element> elements { get; set; }
        public List<BoundaryConditions> boundaryConditions { get; set; }
        
        private double[,] D;
        
        public Elements(Mesh mesh, double[,] D)
        {
            this.D = D;
            elements = new List<Element>();
            boundaryConditions = new List<BoundaryConditions>();
            BuildElements(mesh);
        }

        public class CoordNode //класс для хранения координат
        {
            public CoordNode (double X, double Y, double Z)
            {
                this.X = X;
                this.Y = Y;
                this.Z = Z;
            }
            public double X { get; set; }
            public double Y { get; set; }
            public double Z { get; set; }
        }
        public class Element //класс описывающий элемент
        {
            public Element()
            {
                Node_global = new List<int>();
                coord_Nodes_In_Elemet = new List<CoordNode>();
            }
            public List<int> Node_global { get; set; }
            public List<CoordNode> coord_Nodes_In_Elemet { get; set; }
            public double hx()
            {
                return coord_Nodes_In_Elemet[1].X - coord_Nodes_In_Elemet[0].X;
            }
            public double hy()
            {
                return coord_Nodes_In_Elemet[2].Y - coord_Nodes_In_Elemet[0].Y;
            }
            public double hz()
            {
                return coord_Nodes_In_Elemet[4].Z - coord_Nodes_In_Elemet[0].Z;
            }
            private double Sign(int i)
            {
                if (i == 0) return -1;
                else return 1;
            }
            private double IntegratingFunction(int i, int j, char c1, char c2, double h, double x, bool flag = false)
            {
                double psi1, psi2;
                if (flag)
                {
                    switch (c1)
                    {
                        case 'x':
                            {
                                if (i == 0)
                                    psi1 = (coord_Nodes_In_Elemet[1].X - x) / h;
                                else
                                    psi1 = (x - coord_Nodes_In_Elemet[0].X) / h;

                                return psi1;
                            }
                        case 'y':
                            {
                                switch (c2)
                                {
                                    case 'x':
                                        {
                                            if (i == 0)
                                                psi1 = (coord_Nodes_In_Elemet[2].Y - x) / h;
                                            else
                                                psi1 = (x - coord_Nodes_In_Elemet[0].Y) / h;

                                            return psi1;
                                        }
                                    case 'z':
                                        {
                                            if (i == 0)
                                                psi1 = (coord_Nodes_In_Elemet[1].Y - x) / h;
                                            else
                                                psi1 = (x - coord_Nodes_In_Elemet[0].Y) / h;

                                            return psi1;

                                        }
                                    default: return 0;
                                }
                            }
                        case 'z':
                            {
                                if (i == 0)
                                    psi1 = (coord_Nodes_In_Elemet[2].Z - x) / h;
                                else
                                    psi1 = (x - coord_Nodes_In_Elemet[0].Z) / h;

                                return psi1;
                            }
                        default: return 0;
                    }
                }

                switch (c1)
                {
                    case 'x':
                        {
                            if (i == 0)
                                psi1 = (coord_Nodes_In_Elemet[1].X - x) / h;
                            else
                                psi1 = (x - coord_Nodes_In_Elemet[0].X) / h;
                            break;
                        }
                    case 'y':
                        {
                            if (i == 0)
                                psi1 = (coord_Nodes_In_Elemet[2].Y - x) / h;
                            else
                                psi1 = (x - coord_Nodes_In_Elemet[0].Y) / h;
                            break;
                        }
                    case 'z':
                        {
                            if (i == 0)
                                psi1 = (coord_Nodes_In_Elemet[4].Z - x) / h;
                            else
                                psi1 = (x - coord_Nodes_In_Elemet[0].Z) / h;
                            break;
                        }
                    default:
                        {
                            psi1 = Sign(i) / h;
                            break;
                        }
                }
                switch (c2)
                {
                    case 'x':
                        {
                            if (j == 0)
                                psi2 = (coord_Nodes_In_Elemet[1].X - x) / h;
                            else
                                psi2 = (x - coord_Nodes_In_Elemet[0].X) / h;
                            break;
                        }
                    case 'y':
                        {
                            if (j == 0)
                                psi2 = (coord_Nodes_In_Elemet[2].Y - x) / h;
                            else
                                psi2 = (x - coord_Nodes_In_Elemet[0].Y) / h;
                            break;
                        }
                    case 'z':
                        {
                            if (j == 0)
                                psi2 = (coord_Nodes_In_Elemet[4].Z - x) / h;
                            else
                                psi2 = (x - coord_Nodes_In_Elemet[0].Z) / h;
                            break;
                        }
                    default:
                        {
                            psi2 = Sign(j) / h;
                            break;
                        }
                }
                return psi1 * psi2;
            }
            private double IntegralSimpson(double a, double b, int i, int j, char c1, char c2, double _h, bool flag = false)
            {
                int n = 1;
                double h = b - a;
                double integral = 0;
                for (int k = 0; k < n; k++)
                {
                    double x1 = a + k * h;
                    double x2 = a + (k + 1) * h;

                    if (!flag)
                        integral += ((x2 - x1) / 6) * (IntegratingFunction(i, j, c1, c2, _h, x1) + 4 * IntegratingFunction(i, j, c1, c2, _h, (x1 + x2) / 2) + IntegratingFunction(i, j, c1, c2, _h, x2));
                    else
                        integral += ((x2 - x1) / 6) * (IntegratingFunction(i, j, c1, c2, _h, x1, true) + 4 * IntegratingFunction(i, j, c1, c2, _h, (x1 + x2) / 2, true) + IntegratingFunction(i, j, c1, c2, _h, x2, true));

                }
                return integral;
            }
            public double Kx(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'd', 'd', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'y', 'y', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'z', 'z', hz);
                
                return integralXX * integralYY * integralZZ;
            }
            public double Ky(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'x', 'x', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'd', 'd', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'z', 'z', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kz(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'x', 'x', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'y', 'y', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'd', 'd', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kxy(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'd', 'x', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'y', 'd', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'z', 'z', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kyx(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'x', 'd', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'd', 'y', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'z', 'z', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kxz(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'd', 'x', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'y', 'y', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'z', 'd', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kzx(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'x', 'd', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'y', 'y', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'd', 'z', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kyz(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'x', 'x', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'd', 'y', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'z', 'd', hz);

                return integralXX * integralYY * integralZZ;
            }
            public double Kzy(int i, int j, double hx, double hy, double hz)
            {
                double integralXX = IntegralSimpson(coord_Nodes_In_Elemet[0].X, coord_Nodes_In_Elemet[1].X, mu(i), mu(j), 'x', 'x', hx);
                double integralYY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[2].Y, nu(i), nu(j), 'y', 'd', hy);
                double integralZZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[4].Z, teta(i), teta(j), 'd', 'z', hz);
                
                return integralXX * integralYY * integralZZ;
            }
            private int mu(int i)
            {
                return i % 2;
            }
            private int nu(int i)
            {
                return (i / 2) % 2;
            }
            private int teta(int i)
            {
                return i / 4;
            }
            public double BoundaryIntegral(int i, int Ngran, double hx, double hy, double hz)
            {
                //if (Ngran == 0 || Ngran == 1)
                {
                    double integralY = IntegralSimpson(coord_Nodes_In_Elemet[0].Y, coord_Nodes_In_Elemet[1].Y, mu(i), 0, 'y', 'z', hy, true);
                    double integralZ = IntegralSimpson(coord_Nodes_In_Elemet[0].Z, coord_Nodes_In_Elemet[2].Z, nu(i), 0, 'z', 'z', hz, true);

                    return integralY * integralZ;
                }
                
                //double integral1 = IntegralSimpson(a1, a2, i, j, c1, c2, h1);
                //double integral2 = IntegralSimpson(b1, b2, i, j, c1, c2, h2);

                //return integral1 * integral2;
            }
        }
        public class BoundaryConditions
        {
            public int type { get; set; }
            public int Ngran {  get; set; }
            public List<double> P {  get; set; }
            public double Px {  get; set; }
            public double Py {  get; set; }
            public double Pz {  get; set; }
            public List<Element> elements2D { get; set; }
            public BoundaryConditions(Mesh mesh, int Ngran, int type, double Px, double Py, double Pz)
            {
                elements2D = new List<Element>();
                this.Ngran = Ngran;
                this.type = type;
                this.Px = Px;
                this.Py = Py;
                this.Pz = Pz;

                BuildElements2D(mesh);
            }
            private int GetGlobalNode(int [] paramVec, int NX, int NY)
            {
                return paramVec[0] + (paramVec[1] + paramVec[2] * NY) * NX;
            }
            private void UpdateIndexVec(ref int[] indexVec, bool[] paramVec, int constIndex, int j , int k, bool isFirst, int localIndex)
            {
                if (!isFirst)
                {
                    bool flag = false;
                    for (int i = 0; i < 3; i++)
                        if (paramVec[i])
                        {
                            if (!flag)
                            {
                                flag = true;
                                indexVec[i] = j;
                            }
                            else
                                indexVec[i] = k;
                        }
                        else
                            indexVec[i] = constIndex;
                }
                else
                {
                    bool flag = false;
                    for (int i = 0; i < 3; i++)
                        if (paramVec[i])
                        {
                            if (!flag)
                            {
                                flag = true;
                                indexVec[i] += localIndex % 2;
                            }
                            else
                                indexVec[i] += localIndex / 2;
                        }
                        else
                            indexVec[i] = constIndex;

                }


            }
            private void BuildBoundaryElements(int Nup, int Ndown, int Nx, int Ny, bool[] paramVec, int constIndex, Mesh mesh)
            {
                int[] indexVec = new int[3];
                int[] bufIndexVec = new int[3];

                for (int i = 0; i < 3; i++)
                    if (!paramVec[i])
                        indexVec[i] = constIndex;

                int index, node;
                for (int k = 0; k < Nup - 1; k++)
                    for (int j = 0; j < Ndown - 1; j++)
                    {
                        UpdateIndexVec(ref indexVec, paramVec, constIndex, j, k, false, 0);
                        for (int i = 0; i < 3; i++)
                            bufIndexVec[i] = indexVec[i];

                        index = j + k * (Ndown - 1);

                        elements2D.Add(new Element());

                        for (int localIelem = 0; localIelem < 4; localIelem++)
                        {
                            for (int i = 0; i < 3; i++)
                                indexVec[i] = bufIndexVec[i];

                            UpdateIndexVec(ref indexVec, paramVec, constIndex, j, k, true, localIelem);
                            
                            node = GetGlobalNode(indexVec, Nx, Ny);

                            for (int localIblock = 0; localIblock < 3; localIblock++)
                                elements2D[index].Node_global.Add(3 * node + localIblock);

                            elements2D[index].coord_Nodes_In_Elemet.Add(new CoordNode(mesh.BigMesh[0].isMesh[indexVec[0]], mesh.BigMesh[1].isMesh[indexVec[1]], mesh.BigMesh[2].isMesh[indexVec[2]]));

                        }

                    }

            }
            private void BuildElements2D(Mesh mesh)
            {
                bool[] paramVec;
                switch (Ngran)
                {
                    case 0:
                        {
                            paramVec = new bool[3] { false, true, true};
                            BuildBoundaryElements(mesh.BigMesh[2].isMesh.Count, mesh.BigMesh[1].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[1].isMesh.Count, paramVec, 0, mesh);
                            break;
                        }                    
                    case 1:
                        {
                            paramVec = new bool[3] { false, true, true };
                            BuildBoundaryElements(mesh.BigMesh[2].isMesh.Count, mesh.BigMesh[1].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[1].isMesh.Count, paramVec, mesh.BigMesh[0].isMesh.Count - 1, mesh);
                            break;
                        }
                    case 2:
                        {
                            paramVec = new bool[3] { true, false, true };
                            BuildBoundaryElements(mesh.BigMesh[2].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[1].isMesh.Count, paramVec, 0, mesh);
                            break;
                        }
                    case 3:
                        {
                            paramVec = new bool[3] { true, false, true };
                            BuildBoundaryElements(mesh.BigMesh[2].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[1].isMesh.Count, paramVec, mesh.BigMesh[1].isMesh.Count - 1, mesh);
                            break;
                        }
                    case 4:
                        {
                            paramVec = new bool[3] { true, true, false };
                            BuildBoundaryElements(mesh.BigMesh[1].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[1].isMesh.Count, paramVec, 0, mesh);

                            break;
                        }
                    case 5:
                        {
                            paramVec = new bool[3] { true, true, false };
                            BuildBoundaryElements(mesh.BigMesh[1].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[0].isMesh.Count, mesh.BigMesh[1].isMesh.Count, paramVec, mesh.BigMesh[2].isMesh.Count - 1, mesh);
                            break;
                        }
                    default:
                        {
                            MessageBox.Show("Error! Number of Boundary is not correct!");
                            break;
                        }
                }
            }
        }

        public void BuildBoundaryConditions(Mesh mesh, ObservableCollection<BoundaryCoefs> boundaryCoefs)
        {
            foreach (BoundaryCoefs bc in boundaryCoefs)
            {
                if (bc.Type == 2 && bc.Px == 0 && bc.Py == 0 && bc.Pz == 0)
                    continue;
                
                boundaryConditions.Add(new BoundaryConditions(mesh, bc.Ngran, bc.Type, bc.Px, bc.Py, bc.Pz));
            }
        }
        private double[,] CalcAbloc(int i, int j, Element element)
        {
            double[,] Ablock = new double[3, 3];

            double hx = element.hx(),
                hy = element.hz(),
                hz = element.hz(),
                Kx = element.Kx(i, j, hx, hy, hz),
                Ky = element.Ky(i, j, hx, hy, hz),
                Kz = element.Kz(i, j, hx, hy, hz),
                Kxy = element.Kxy(i, j, hx, hy, hz),
                Kxz = element.Kxz(i, j, hx, hy, hz),
                Kyx = element.Kyz(i, j, hx, hy, hz),
                Kyz = element.Kyz(i, j, hx, hy, hz),
                Kzx = element.Kzx(i, j, hx, hy, hz),
                Kzy = element.Kzy(i, j, hx, hy, hz);

            Ablock[0, 0] = D[0, 0] * Kx + D[3, 3] * Ky + D[5, 5] * Kz;
            Ablock[0, 1] = D[0, 1] * Kxy + D[3, 3] * Kyx;
            Ablock[0, 2] = D[0, 2] * Kxz + D[3, 3] * Kzx;

            Ablock[1, 0] = D[1, 0] * Kyx + D[3, 3] * Kxy;
            Ablock[1, 1] = D[3, 3] * Kx + D[1, 1] * Ky + D[4, 4] * Kz;
            Ablock[1, 2] = D[1, 2] * Kyz + D[4, 4] * Kzy;

            Ablock[2, 0] = D[2, 0] * Kzx + D[5, 5] * Kxz;
            Ablock[2, 1] = D[2, 1] * Kzy + D[4, 4] * Kyz;
            Ablock[2, 2] = D[5, 5] * Kx + D[4, 4] * Ky + D[2,2] * Kz;

            return Ablock;
        }
        public void BuildElements(Mesh mesh)
        {
            int Nx = mesh.BigMesh[0].isMesh.Count - 1,
                Ny = mesh.BigMesh[1].isMesh.Count - 1,
                Nz = mesh.BigMesh[2].isMesh.Count - 1,
                index, node;
            
            for (int k = 0; k < Nz; k++)
                for (int j = 0; j < Ny; j++)
                    for (int i = 0; i < Nx; i++)
                    {
                        index = i + (j + k * Ny) * Nx;
                        elements.Add(new Element());

                        for (int localIelem = 0; localIelem < 8; localIelem++)
                        {
                            for (int localIblock  = 0; localIblock < 3;  localIblock++)
                            {
                                node = 3 * ((i + localIelem % 2) + ((j + (localIelem / 2) % 2) + (k + localIelem / 4) * (Ny + 1)) * (Nx + 1)) + localIblock;
                                elements[index].Node_global.Add(node);
                            }
                            elements[index].coord_Nodes_In_Elemet.Add(new CoordNode(
                                mesh.BigMesh[0].isMesh[i + localIelem % 2],
                                mesh.BigMesh[1].isMesh[j + (localIelem / 2) % 2],
                                mesh.BigMesh[2].isMesh[k + localIelem / 4]));
                        }
                    }
        }
        public void BuildGlobalMatrix(int[] ig, int[] jg, ref double[] di, ref double[] ggl)
        {
            double[,] Abloc;

            foreach (Element element in elements)
            {
                List<List<double>> A_loc = new List<List<double>>();
                for (int i = 0; i < 24; i++)
                {
                    A_loc.Add(new List<double>());
                    for (int j = 0; j < 24; j++)
                        A_loc[i].Add(0);
                }

                for (int i = 0; i < 8; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        Abloc = CalcAbloc(i, j, element);

                        for (int ii = 0; ii < 3; ii++)
                            for (int jj = 0; jj < 3; jj++)
                                A_loc[3 * i + ii][3 * j + jj] += Abloc[ii, jj];
                    }
                }
                AddToMatrix(A_loc, element, ig, jg, ref di, ref ggl);
            }
        }
        private void AddToMatrix(List<List<double>> A_loc, Element element, int[] ig, int[] jg, ref double[] di, ref double[] ggl)//функция для внесения локальных элементов матрицы в глобальную
        {
            List<int> L = element.Node_global;
            int n_loc = element.Node_global.Count;

            for (int i = 0; i < n_loc; i++)
            {
                di[L[i]] += A_loc[i][i];
            }

            for (int i = 0; i < n_loc; i++)
            {
                int temp = ig[L[i]];
                for (int j = 0; j < i; j++)
                    for (int k = temp; k < ig[L[i] + 1]; k++)
                    {
                        if (jg[k] == L[j])
                        {
                            ggl[k] += A_loc[i][j];
                            break;
                        }
                    }
            }
        }

    }
    public class Mesh
    {
        public class OneDimensionalMesh
        {
            public List<double> isMesh { get; set; }
            public OneDimensionalMesh(double h, double kRazr, double coordEnd)
            {
                isMesh = new List<double>();
                BuildMesh(h, kRazr, coordEnd);
            }
            private void BuildMesh(double h, double kRazr, double coordEnd)
            {
                isMesh.Clear();
                int i = 0;
                isMesh.Add(0);

                if (h >= coordEnd)
                    isMesh.Add(coordEnd);
                else
                {
                    while (isMesh[i] < coordEnd)
                    {
                        i++;
                        isMesh.Add(isMesh[i - 1] + h * Math.Pow(kRazr, i));
                    }

                    if ((isMesh[i - 1] + isMesh[i]) / 2 > coordEnd)
                    {
                        isMesh[i - 1] = coordEnd;
                        isMesh.RemoveAt(i);
                    }
                    else
                        isMesh[i] = coordEnd;
                }
            }
        }
        public List<OneDimensionalMesh> BigMesh { get; set; }
        private double SizeX, SizeY, SizeZ;
        private MeshCoefs Coefs;
        public Mesh(double SizeX, double SizeY, double SizeZ, MeshCoefs Coefs) 
        {
            this.Coefs = Coefs;
            this.SizeX = SizeX;
            this.SizeY = SizeY;
            this.SizeZ = SizeZ;
            BigMesh = new List<OneDimensionalMesh>();
        }
        public void BuildMesh3D()
        {
            BigMesh.Clear();

            BigMesh.Add(new OneDimensionalMesh(Coefs.hX, Coefs.kX, SizeX));
            BigMesh.Add(new OneDimensionalMesh(Coefs.hY, Coefs.kY, SizeY));
            BigMesh.Add(new OneDimensionalMesh(Coefs.hZ, Coefs.kZ, SizeZ));
        }
    }
}
