namespace Actuators
{
    /// <summary>
    /// Класс ЭМ
    /// </summary>
    public class Magnetorquer
    {
        /* включение ЭМ:
         * 1 - включен в положительном направлении
         * -1 - включен в отрицательном направлении
         * 0 - выключен
        */
        private int activation = 0;
        private double[] L = new double[3];  // [А*м^2] - магнитный момент ЭМ


        /// <summary>
        /// Название ЭМ
        /// </summary>
        public string Name { get; set; }

        /// <summary>
        /// Магнитный момент ЭМ, [А*м^2]
        /// </summary>
        public double MagneticMomentNominal { get; set; }

        /// <summary>
        /// Орт установки ЭМ относительно ССК (соответствует положительному магнитному моменту)
        /// </summary>
        public double[] Orient { get; set; } = new double[3];

        /// <summary>
        /// Функционирование ЭМ (true - работает, false - не работает)
        /// </summary>
        public bool OperationalStatus { get; set; }

        /// <summary>
        /// Флаг включения (+1, -1 - включен, 0 - выключен)
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
                // магнитный момент
                L = MagneticMomentCalc(activation);
            }
            get { return activation; }
        }

        /// <summary>
        /// Расчет проекций магнитного момента ЭМ на оси ССК
        /// </summary>
        /// <param name="activation">Включение ЭМ (+1, -1 - включен, 0 - выключен)</param>
        /// <returns>Проекции создаваемого магнитного момента на оси ССК</returns>
        private double[] MagneticMomentCalc(int activation)
        {
            if (activation == 0)
                return new double[3] { 0, 0, 0 };
            else
                return new double[3] { activation * MagneticMomentNominal * Orient[0],
                                       activation * MagneticMomentNominal * Orient[1],
                                       activation * MagneticMomentNominal * Orient[2] };
        }

        /// <summary>
        ///  Проекции магнитного момента ЭМ на оси ССК, [А*м^2]
        /// </summary>
        public double[] GetMagneticMoment
        {
            get { return L; }
        }

        /// <summary>
        ///  Проекции магнитного момента ЭМ на оси ССК, [А*м^2]
        /// </summary>
        public double[] MagneticMoment
        {
            get { return L; }
        }

        /// <summary>
        /// Момент, создаваемый ЭМ
        /// </summary>
        /// <param name="B">Магнитная индукция МПЗ в ССК, [Тл]</param>
        /// <returns>Проеции создаваемого момента на оси ССК</returns>
        public double[] GetTorque(double[] B)
        {
            return new double[3] { L[1] * B[2] - L[2] * B[1],
                                   L[2] * B[0] - L[0] * B[2],
                                   L[0] * B[1] - L[1] * B[0] };
        }

        /// <summary>
        /// Момент, создаваемый ЭМ
        /// </summary>
        /// <param name="B">Магнитная индукция МПЗ в ССК, [Тл]</param>
        /// <returns>Проеции создаваемого момента на оси ССК</returns>
        public double[] Torque(double[] B)
        {
            return new double[3] { L[1] * B[2] - L[2] * B[1],
                                   L[2] * B[0] - L[0] * B[2],
                                   L[0] * B[1] - L[1] * B[0] };
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="name">Название ЭМ</param>
        /// <param name="magneticMomentNominal">Магнитный момент ЭМ, [А*м^2]</param>
        /// <param name="orientX">X-компонента орта установки ЭМ относительно ССК</param>
        /// <param name="orientY">Y-компонента орта установки ЭМ относительно ССК</param>
        /// <param name="orientZ">Z-компонента орта установки ЭМ относительно ССК</param>
        /// <param name="operationalStatus">Функционирование ЭМ (true - работает, false - не работает)</param>
        public Magnetorquer(string name = "ЭМ",
                            double magneticMomentNominal = 1,
                            double orientX = 1, double orientY = 0, double orientZ = 0,
                            bool operationalStatus = true)
        {
            // название ЭМ
            Name = name;
            // номинальный магнитный момент
            MagneticMomentNominal = magneticMomentNominal;
            // ориентация ЭМ
            Orient[0] = orientX;
            Orient[1] = orientY;
            Orient[2] = orientZ;
            // функционирование ЭМ
            OperationalStatus = operationalStatus;
        }
    }
}
