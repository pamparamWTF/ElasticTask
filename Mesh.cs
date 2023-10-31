using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace ElasticTask
{
    public class SLAE
    {

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
        //private GeometryData geometryData;
        private double[,] D;
        public Elements(Mesh mesh, double[,] D)
        {
            this.D = D;
            elements = new List<Element>();
            BuildElements(mesh);
        }
        public class Coord_Node //класс для хранения координат
        {
            public Coord_Node (double X, double Y, double Z)
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
                coord_Nodes_In_Elemet = new List<Coord_Node>();
            }
            public List<int> Node_global { get; set; }
            public List<Coord_Node> coord_Nodes_In_Elemet { get; set; }
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
            private double IntegratingFunction(int i, int j, char c1, char c2, double h, double x)
            {
                double psi1, psi2;
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
            private double IntegralSimpson(double a, double b, int i, int j, char c1, char c2, double _h)
            {
                int n = 1;
                double h = b - a;
                double integral = 0;
                for (int k = 0; k < n; k++)
                {
                    double x1 = a + k * h;
                    double x2 = a + (k + 1) * h;

                    integral += ((x2 - x1) / 6) * (IntegratingFunction(i, j, c1, c2, _h, x1) + 4 * IntegratingFunction(i, j, c1, c2, _h, (x1 + x2) / 2) + IntegratingFunction(i, j, c1, c2, _h, x2));
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
                            elements[index].coord_Nodes_In_Elemet.Add(new Coord_Node(
                                mesh.BigMesh[0].isMesh[i + localIelem % 2],
                                mesh.BigMesh[1].isMesh[j + (localIelem / 2) % 2],
                                mesh.BigMesh[2].isMesh[k + localIelem / 4]));
                        }
                    }
        }
        public void BuildGlobalMatrix()
        {
            int n = elements.Count();
            double[,] Abloc = new double[3,3];

            for (int elem = 0; elem < n; elem++)
            {
                List<List<double>> A_loc = new List<List<double>>();
                for (int i = 0; i < 8; i++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        A_loc.Add(new List<double>());
                    
                    for (int j = 0; j < 8; j++)
                    {
                        Abloc = CalcAbloc(i, j, elements[elem]);
                        
                        for (int ii = 0; ii < 3; ii++)
                            for (int jj = 0; jj < 3; jj++) 
                                A_loc[i + ii].Add(Abloc[ii, jj]);
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
