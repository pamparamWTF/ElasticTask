using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
using System.Windows.Shapes;

namespace ElasticTask
{
    /// <summary>
    /// Логика взаимодействия для GeometryWindow.xaml
    /// </summary>
    /// 
    public partial class GeometryWindow : Window
    {
        public ObservableCollection<GeometryData> geometryData {  get; private set; }
        public GeometryWindow(ref ObservableCollection<GeometryData> geometryData)
        {
            InitializeComponent();
            this.SizeToContent = SizeToContent.Width;
            this.geometryData = geometryData;
            this.DataContext = geometryData;
            GeometryData.ItemsSource = this.geometryData;

        }
        private void AddGemetry_Click(object sender, RoutedEventArgs e)
        {
            geometryData.Add(new GeometryData(geometryData.Count + 1));
            GeometryData.Items.Refresh();
        }
        private void DeleteGeometry_Click(object sender, RoutedEventArgs e)
        {
            if (GeometryData.SelectedIndex != -1)
            {
                geometryData.RemoveAt(GeometryData.SelectedIndex);
                for (int i = 1; i <= geometryData.Count; i++) geometryData[i - 1].N = i;
                GeometryData.Items.Refresh();
            }
        }
        private void CopyGeometry_Click(object sender, RoutedEventArgs e)
        {
            GeometryData copied = new GeometryData(GeometryData.SelectedIndex)
            {
                NameGeometry = geometryData[GeometryData.SelectedIndex].NameGeometry,
                SizeX = geometryData[GeometryData.SelectedIndex].SizeX, SizeY = geometryData[GeometryData.SelectedIndex].SizeY, SizeZ = geometryData[GeometryData.SelectedIndex].SizeZ,
                NameMaterial = geometryData[GeometryData.SelectedIndex].NameMaterial,
                Ex = geometryData[GeometryData.SelectedIndex].Ex,
                Ey = geometryData[GeometryData.SelectedIndex].Ey,
                Ez = geometryData[GeometryData.SelectedIndex].Ez,
                Vxy = geometryData[GeometryData.SelectedIndex].Vxy,
                Vxz = geometryData[GeometryData.SelectedIndex].Vxz,
                Vyz = geometryData[GeometryData.SelectedIndex].Vyz,
                Gxy = geometryData[GeometryData.SelectedIndex].Gxy,
                Gxz = geometryData[GeometryData.SelectedIndex].Gxz,
                Gyz = geometryData[GeometryData.SelectedIndex].Gyz,
                OnModel = false
            };
            geometryData.Add(copied);
            for (int i = 1; i <= geometryData.Count; i++) geometryData[i - 1].N = i;
            GeometryData.Items.Refresh();
        }
    }
}