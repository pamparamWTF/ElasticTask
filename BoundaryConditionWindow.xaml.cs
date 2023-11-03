using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
    /// Логика взаимодействия для BoundaryConditionWindow.xaml
    /// </summary>
    public partial class BoundaryConditionWindow : Window
    {
        public ObservableCollection<BoundaryCoefs> ListOfBoundaryCoefs { get; private set; }
        public BoundaryConditionWindow(ref ObservableCollection<BoundaryCoefs> ListOfBoundaryCoefs)
        {
            InitializeComponent();

            this.ListOfBoundaryCoefs = ListOfBoundaryCoefs;
            this.SizeToContent = SizeToContent.Width;
            this.DataContext = ListOfBoundaryCoefs;
            BoundaryData.ItemsSource = this.ListOfBoundaryCoefs;
        }
    }
}
