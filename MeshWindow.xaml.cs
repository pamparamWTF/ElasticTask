using System;
using System.Collections.Generic;
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
    /// Логика взаимодействия для MeshWindow.xaml
    /// </summary>
    public partial class MeshWindow : Window
    {
        public MeshCoefs meshCoefs { get; private set; }
        public MeshWindow(ref MeshCoefs meshCoefs)
        {
            InitializeComponent();
            this.meshCoefs = meshCoefs;
            this.DataContext = this.meshCoefs;
        }
    }
}
