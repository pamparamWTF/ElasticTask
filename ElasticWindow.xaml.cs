using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.IO;
using System.Windows.Forms;
using Newtonsoft.Json;
using Microsoft.Win32;
using System.Collections.ObjectModel;

namespace ElasticTask
{
    /// <summary>
    /// Interaction logic for ElasticWindow.xaml
    /// </summary>
    public partial class ElasticWindow : Window
    {
        string path = "";
        GeometryWindow GeometryWindow;
        MeshWindow MeshWindow;
        Elements elements;
        public struct Model
        {
            public Model(ObservableCollection<GeometryData> ListOfGeometry, Mesh mesh, MeshCoefs meshCoefs)
            {
                this.ListOfGeometry = ListOfGeometry;
                this.mesh = mesh;
                this.meshCoefs = meshCoefs;
            }
            public ObservableCollection<GeometryData> ListOfGeometry;
            public Mesh mesh;
            public MeshCoefs meshCoefs;
        }
        Model model = new(new ObservableCollection<GeometryData>() { new GeometryData(1) }, new Mesh(1, 1, 1, new MeshCoefs(1, 1, 1, 1, 1, 1)), new MeshCoefs(1, 1, 1, 1, 1, 1));

        public ElasticWindow()
        {
            InitializeComponent();
            model.ListOfGeometry[0].OnModel = true;
        }
        private void ClickMenuMesh(object sender, RoutedEventArgs e)
        {
            MeshWindow = new MeshWindow(ref model.meshCoefs);
            bool? diaRes = MeshWindow.ShowDialog();

            if (diaRes.HasValue && diaRes.Value)
            {
                model.meshCoefs = MeshWindow.meshCoefs;
            }

        }
        private void ClickMenuGeometry(object sender, RoutedEventArgs e)
        {
            GeometryWindow = new GeometryWindow(ref model.ListOfGeometry);
            bool? diaRes = GeometryWindow.ShowDialog();

            if (diaRes.HasValue && diaRes.Value)
            {
                model.ListOfGeometry = GeometryWindow.geometryData;
            }
        }
        private void SaveModelAs()
        {
            System.Windows.Forms.SaveFileDialog saveFileDialog = new System.Windows.Forms.SaveFileDialog();
            saveFileDialog.Filter = "Elastic Task Model files (*.etm)|*.etm";
            
            if (saveFileDialog.ShowDialog() == System.Windows.Forms.DialogResult.Cancel)
                return;

            path = saveFileDialog.FileName;

            if (path == "")
                System.Windows.Forms.MessageBox.Show("Error! Could not save model.", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            else
            {
                SaveModel(path); 
            }
        }
        private void SaveModel(string path)
        {
            var jsonString = JsonConvert.SerializeObject(model);
            File.WriteAllText(path, jsonString);
        }
        private void ClickMenuOpen(object sender, RoutedEventArgs e)
        {
            System.Windows.Forms.OpenFileDialog openFileDialog = new System.Windows.Forms.OpenFileDialog();
            openFileDialog.Filter = "Elastic Task Model files (*.etm)|*.etm";
            
            if (openFileDialog.ShowDialog() == System.Windows.Forms.DialogResult.Cancel)
                return;

            path = openFileDialog.FileName;

            model = JsonConvert.DeserializeObject<Model>(File.ReadAllText(path));
        }
        private void ClickMenuSave(object sender, RoutedEventArgs e)
        {
            if (path == "")
                SaveModelAs();
            else
                SaveModel(path);
        }
        private void ClickMenuSaveAs(object sender, RoutedEventArgs e) => SaveModelAs();
        private void ClickMenuExit(object sender, RoutedEventArgs e) => Close();
        protected override void OnClosing(CancelEventArgs e)
        {
            if (System.Windows.MessageBox.Show("Are you sure you want to quit?",
                    "Save file",
                    MessageBoxButton.YesNo,
                    MessageBoxImage.Question) == MessageBoxResult.Yes)
                System.Windows.Application.Current.Shutdown();
            else e.Cancel = true;
        }
        private void BuildMeshButton_Click(object sender, RoutedEventArgs e)
        {
            if (path == "")
                SaveModelAs();

            if (path == "")
                return;

            int count = 0, N = 0;
            foreach (var geometry in model.ListOfGeometry)
                if (geometry.OnModel) { count++; N = geometry.N - 1; }

            if (count != 1)
                System.Windows.Forms.MessageBox.Show("Error! Could not build mesh. Please check \"Model -> Geometry\" and enter only one geometry.", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            else
            {
                model.mesh = new Mesh(
                    model.ListOfGeometry[N].SizeX,
                    model.ListOfGeometry[N].SizeY,
                    model.ListOfGeometry[N].SizeZ,
                    model.meshCoefs);

                model.mesh.BuildMesh3D();

                elements = new Elements(model.mesh);

                //System.Windows.Forms.MessageBox.Show(elements.elements[0].Ky(0, 0, elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                ////+ " " +
                ////    elements.elements[0].Kx(0, 1).ToString() + " " +
                ////    elements.elements[0].Kx(0, 2).ToString() + " " +
                ////    elements.elements[0].Kx(0, 3).ToString() + " " +
                ////    elements.elements[0].Kx(0, 4).ToString() + " " +
                ////    elements.elements[0].Kx(0, 5).ToString() + " " +
                ////    elements.elements[0].Kx(0, 6).ToString() + " " +
                ////    elements.elements[0].Kx(0, 7).ToString());

                SaveModel(path);
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            string s = tb1.Text + tb2.Text;
            switch (s)
            {
                case "xx":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kx(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "yy":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Ky(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "zz":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kz(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "xy":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kxy(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "xz":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kxz(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "yx":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kyx(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "yz":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kyz(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "zx":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kzx(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
                case "zy":
                    {
                        System.Windows.Forms.MessageBox.Show(elements.elements[0].Kzy(int.Parse(num1.Text), int.Parse(num2.Text), elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
                        break;
                    }
            }
            //System.Windows.Forms.MessageBox.Show(elements.elements[0].Ky(0, 0, elements.elements[0].hx(), elements.elements[0].hy(), elements.elements[0].hz()).ToString());
            //+ " " +
            //    elements.elements[0].Kx(0, 1).ToString() + " " +
            //    elements.elements[0].Kx(0, 2).ToString() + " " +
            //    elements.elements[0].Kx(0, 3).ToString() + " " +
            //    elements.elements[0].Kx(0, 4).ToString() + " " +
            //    elements.elements[0].Kx(0, 5).ToString() + " " +
            //    elements.elements[0].Kx(0, 6).ToString() + " " +
            //    elements.elements[0].Kx(0, 7).ToString());

        }
    }
}
