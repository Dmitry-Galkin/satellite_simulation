namespace Actuators
{
    /// <summary>
    /// Класс ЖРД
    /// </summary>
    public class Thruster
    {
        private int activation = 0;  // открытие клапана (1 - открыт, 0 - закрыт)
        private double[] thrust = new double[3];  // [Н], проекции тяги на оси ССК
        private double[] torque = new double[3];  // [Н*м], проекции механического момента на оси ССК


        /// <summary>
        /// Название двигателя
        /// </summary>
        public string Name { get; set; }

        /// <summary>
        /// Тяга двигателя, [Н]
        /// </summary>
        public double ThrustNominal { get; set; }

        /// <summary>
        /// Координаты установки двигателя относительно ЦМ в ССК, [м]
        /// </summary>
        public double[] Coord { get; set; } = new double[3];

        /// <summary>
        /// Орт ориентации вектора тяги двигателя относительно ССК
        /// </summary>
        public double[] Orient { get; set; } = new double[3];

        /// <summary>
        /// Функционирование двигателя (true - работает, false - не работает)
        /// </summary>
        public bool OperationalStatus { get; set; }

        /// <summary>
        /// Флаг открытия клапана (1 - открыт, 0 - закрыт)
        /// </summary>
        public int Activation
        {
            set
            {
                if (OperationalStatus)
                {
                    activation = value;
                }
                else
                {
                    activation = 0;
                }
                // тяга
                thrust = ThrustCalc(activation);
                // момент
                torque = TorqueCalc(thrust);
            }
            get { return activation; }
        }

        /// <summary>
        /// Расчет проекций создаваемой двигателем тяги на оси ССК
        /// </summary>
        /// <param name="activation">Флаг открытия клапана (1 - открыт, 0 - закрыт)</param>
        /// <returns>Проекции создаваемой тяги на оси ССК</returns>
        private double[] ThrustCalc(int activation)
        {
            if (activation == 0)
                return new double[3] { 0, 0, 0 };
            else
                return new double[3] { ThrustNominal * Orient[0], 
                                       ThrustNominal * Orient[1], 
                                       ThrustNominal * Orient[2] };
        }

        /// <summary>
        /// Расчет проекций создаваемого двигателем момента на оси ССК
        /// </summary>
        /// <param name="thrust">Вектор создаваемой двигателем тяги</param>
        /// <returns>Проекции создаваемого момента на оси ССК</returns>
        private double[] TorqueCalc(double[] thrust)
        {
            return new double[3] { Coord[1] * thrust[2] - Coord[2] * thrust[1],
                                   Coord[2] * thrust[0] - Coord[0] * thrust[2],
                                   Coord[0] * thrust[1] - Coord[1] * thrust[0] };
        }

        /// <summary>
        /// Проекции создаваемой двигателем тяги на оси ССК, [Н]
        /// </summary>
        public double[] GetThrust
        { 
            get { return thrust; } 
        }

        /// <summary>
        /// Проекции создаваемой двигателем тяги на оси ССК, [Н]
        /// </summary>
        public double[] Thrust
        {
            get { return thrust; }
        }

        /// <summary>
        /// Проекции создаваемого двигателем момента на оси ССК, [Н*м]
        /// </summary>
        public double[] GetTorque
        {
            get { return torque; }
        }

        /// <summary>
        /// Проекции создаваемого двигателем момента на оси ССК, [Н*м]
        /// </summary>
        public double[] Torque
        {
            get { return torque; }
        }

        /// <summary>
        /// Конструктор класса ЖРД
        /// </summary>
        /// <param name="name">Название двигателя</param>
        /// <param name="thrustNominal">Тяга двигателя, [Н]</param>
        /// <param name="coordX">Координата X установки двигателя относительно ЦМ в ССК, [м]</param>
        /// <param name="coordY">Координата Y установки двигателя относительно ЦМ в ССК, [м]</param>
        /// <param name="coordZ">Координата Z установки двигателя относительно ЦМ в ССК, [м]</param>
        /// <param name="orientX">X-компонента орта ориентации вектора тяги относительно ССК</param>
        /// <param name="orientY">Y-компонента орта ориентации вектора тяги относительно ССК</param>
        /// <param name="orientZ">Z-компонента орта ориентации вектора тяги относительно ССК</param>
        /// <param name="operationalStatus">Функционирование двигателя (true - работает, false - не работает)</param>
        public Thruster(string name = "ЖРД",
                        double thrustNominal = 1,
                        double coordX = 0, double coordY = 0, double coordZ = 0,
                        double orientX = 1, double orientY = 0, double orientZ = 0,
                        bool operationalStatus = true)
        {
            // название двигателя
            Name = name;
            // номинальная тяга
            ThrustNominal = thrustNominal;
            // координаты установки
            Coord[0] = coordX;
            Coord[1] = coordY;
            Coord[2] = coordZ;
            // ориентация вектора тяги
            Orient[0] = orientX;
            Orient[1] = orientY;
            Orient[2] = orientZ;
            // функционирование двигателя
            OperationalStatus = operationalStatus;
        }
    }
}
