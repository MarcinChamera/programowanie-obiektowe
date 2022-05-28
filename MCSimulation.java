import java.util.ArrayList;
import java.util.List;
import java.lang.Math;
import java.util.Random;

public class MCSimulation implements main.Simulation {

    private LatticeParametersImpl latticeParametersImpl = new LatticeParametersImpl();
    private ProbabilityFormula formula;
    private double TkB;
    private double Ce;
    private ArrayList<Double> Cn;
    private double externalFieldAngle;
    private int magnetsCount;

    public MCSimulation() {

    }

    @Override
    public void setLattice(int[][] lattice, int states) {
        latticeParametersImpl.setStates(states);
        latticeParametersImpl.setLattice(lattice);
        magnetsCount = lattice.length * lattice.length;
    }

    @Override
    public void setEnergyParameters(List<Double> parameters, double externaFieldAngle) {
        Ce = parameters.get(0);
        // Tutaj dla Cn element z zerowym indeksem to Ce. C1, C2 itd są kolejno z pierwszym, drugim, ... indeksem
        Cn = (ArrayList<Double>) parameters;
        externalFieldAngle = externaFieldAngle;
    }

    @Override
    public void setProbabilityFormula(ProbabilityFormula formula) {
        this.formula = formula;
    }

    @Override
    public void setTkB(double TkB) {
        this.TkB = TkB;
    }

    @Override
    public LatticeParameters getState() {
        return latticeParametersImpl;
    }

    @Override
    public void executeMCSteps(int steps) {
        double totalEnergy = calculateTotalEnergy(latticeParametersImpl.lattice());
        latticeParametersImpl.setTotalEnergy(totalEnergy);
        double acceptances = 0;
        double acceptanceRatio = 0;
        // double orderParameter = 0;
        Random random; 
        for (int step = 0; step < steps; step++) {
            // orderParameter += calculateOrderParameter();
            // System.out.println(orderParameter); 
            if (step > 0) {
                acceptanceRatio = acceptances / step;
            }
            // System.out.println(acceptanceRatio); 
            double deltaE = 0;
            random = new Random();
            int magnetRowRandom = random.nextInt((int)Math.sqrt(magnetsCount));
            int magnetColRandom = random.nextInt((int)Math.sqrt(magnetsCount));
            int[][] newLattice = generateLatticeCopy(latticeParametersImpl.lattice());
            if (acceptanceRatio > 0.5) {
                boolean moreThanOneMagnet = random.nextBoolean();
                // Jeśli są tylko dwa dostępne kąty dla magnesu, to nie ma sensu zmiana tego kąta o więcej niż jeden - dlatego należy
                // zmienić dwa magnesy.
                if (latticeParametersImpl.states() == 2) {
                    moreThanOneMagnet = true;
                }
                if (moreThanOneMagnet) {
                    int magnetRowRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                    int magnetColRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                    while (magnetRowRandom == magnetRowRandom2 && magnetColRandom == magnetColRandom2) {
                        magnetRowRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                        magnetColRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                    }
                    int magnetStateChange = random.nextBoolean() ? 1 : -1;
                    int magnetStateChange2 = random.nextBoolean() ? 1 : -1;
                    changeMagnetState(newLattice, magnetRowRandom, magnetColRandom, magnetStateChange);
                    changeMagnetState(newLattice, magnetRowRandom2, magnetColRandom2, magnetStateChange2);
                    deltaE = calculateTotalEnergy(newLattice) - calculateTotalEnergy(latticeParametersImpl.lattice());
                } else {
                    int magnetStateChange = random.nextBoolean() ? 2 : -2;
                    changeMagnetState(newLattice, magnetRowRandom, magnetColRandom, magnetStateChange);
                    deltaE = calculateEi(newLattice, magnetRowRandom, magnetColRandom);
                }
            } else {
                int magnetStateChange = random.nextBoolean() ? 1 : -1;
                changeMagnetState(newLattice, magnetRowRandom, magnetColRandom, magnetStateChange);
                deltaE = calculateEi(newLattice, magnetRowRandom, magnetColRandom);
            }
            double R = Math.random();
            double P = calculateP(formula, deltaE, TkB);
            if (R < P) {
                totalEnergy += deltaE;
                latticeParametersImpl.setTotalEnergy(totalEnergy);
                latticeParametersImpl.setLattice(newLattice);
                acceptances++;
            }
        }
        // System.out.println("Order parameter srednia: " + orderParameter / steps);
    }

    private double calculateEi(int[][] lattice, int i_row, int i_col) {
        double Ei = 0;
        int tmp_neighboursStatesSize = 0;
        int tmp_j = 0;
        try {
            for (int n = 0; n < magnetsCount; n++) {
                ArrayList<Integer> neighboursStates = getNeighboursStates(lattice, i_row, i_col, n);
                tmp_neighboursStatesSize = neighboursStates.size();
                for (int j = 0; j < neighboursStates.size(); j++) { 
                    tmp_j = j;
                    double alphaI = getAngleInRadians(latticeParametersImpl.magnetState(i_row, i_col));
                    double alphaJ = getAngleInRadians(neighboursStates.get(j));
                    double current_Cn = 0;
                    if (n < Cn.size()) {
                        current_Cn = Cn.get(n);
                    }
                    Ei -= current_Cn * Math.cos(alphaI - alphaJ);
                }
            }
            // odtąd
            Ei *= 0.5;
        } catch(Exception exception) {
            System.out.println("Exception: " + exception + " dla i_row=" + Double.toString(i_row) + " i_col=" + Double.toString(i_col) + " tmp_neighboursStatesSize=" + Integer.toString(tmp_neighboursStatesSize) + " tmp_j=" + Integer.toString(tmp_j));
        }
            
        double alphaI = getAngleInRadians(latticeParametersImpl.magnetState(i_row, i_col));
        Ei -= Ce * Math.cos(alphaI - externalFieldAngle);
        // dotąd - nie wiem czy ma to być czy nie
        return Ei;
    }

    private double calculateTotalEnergy(int[][] lattice) {
        double totalEnergy = -0.5;
        double nSum = 0;
        for (int n = 1; n < Cn.size(); n++) {
            double iSum = 0;
            for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
                for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                    ArrayList<Integer> neighboursStates = getNeighboursStates(lattice, i_row, i_col, n);
                    double jSum = 0;
                    for (int j = 0; j < neighboursStates.size(); j++) {
                        double alphaI = getAngleInRadians(latticeParametersImpl.magnetState(i_row, i_col));
                        double alphaJ = getAngleInRadians(neighboursStates.get(j));
                        jSum +=  Math.cos(alphaI - alphaJ);
                    }
                    iSum += jSum;
                }
            }
            iSum *= Cn.get(n);
            nSum += iSum;
        }
        totalEnergy *= nSum;
        double subtract = Ce;
        double iSum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                double alphaI = getAngleInRadians(latticeParametersImpl.magnetState(i_row, i_col));
                iSum += Math.cos(alphaI - externalFieldAngle);
            }
        }
        subtract *= iSum;
        totalEnergy -= subtract;
        return totalEnergy;
    }

    private ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol, int level) {
        ArrayList<Integer> neighboursStates = new ArrayList<Integer>();
        if (level == 1) {
            if (magnetRow > 0) neighboursStates.add(lattice[magnetRow - 1][magnetCol]);
            if (magnetRow < lattice.length - 1) neighboursStates.add(lattice[magnetRow + 1][magnetCol]);
            if (magnetCol > 0) neighboursStates.add(lattice[magnetRow][magnetCol - 1]);
            if (magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow][magnetCol + 1]);

        } else if (level == 2) {
            if (magnetRow > 0 && magnetCol > 0) neighboursStates.add(lattice[magnetRow - 1][magnetCol - 1]);
            if (magnetRow < lattice.length - 1 && magnetCol > 0) neighboursStates.add(lattice[magnetRow + 1][magnetCol - 1]);
            if (magnetRow > 0 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow][magnetCol - 1]);
            if (magnetRow < lattice.length - 1 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow][magnetCol + 1]);

        } else if (level == 3) {
            if (magnetRow > 1) neighboursStates.add(lattice[magnetRow - 2][magnetCol]);
            if (magnetRow < lattice.length - 2) neighboursStates.add(lattice[magnetRow + 2][magnetCol]);
            if (magnetCol > 1) neighboursStates.add(lattice[magnetRow][magnetCol - 2]);
            if (magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow][magnetCol + 2]);

        } else if (level == 4) {
            if (magnetRow > 0 && magnetCol > 1) neighboursStates.add(lattice[magnetRow - 1][magnetCol - 2]);
            if (magnetRow < lattice.length - 1 && magnetCol > 1) neighboursStates.add(lattice[magnetRow + 1][magnetCol - 2]);
            if (magnetRow > 0 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow][magnetCol - 2]);
            if (magnetRow < lattice.length - 1 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow][magnetCol + 2]);
        }
        // TODO jak to zrobić dla wszystkich poziomów?
        return neighboursStates;
    }

    private double getAngleInRadians(int magnetState) {
        return 2 * Math.PI * magnetState / latticeParametersImpl.states();
    }

    private int[][] generateLatticeCopy(int[][] lattice) {
        int[][] latticeCopy = new int[lattice.length][lattice.length];
        for (int i = 0; i < lattice.length; i++) {
            for (int j = 0; j < lattice.length; j++) {
                latticeCopy[i][j] = lattice[i][j];
            }
        }
        return latticeCopy;
    }

    private double calculateOrderParameter() {
        int[][] lattice = latticeParametersImpl.lattice();
        double xAvg = 1. / magnetsCount;
        double sum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                sum += Math.cos(getAngleInRadians(lattice[i_row][i_col]));
            }
        }
        xAvg *= sum;
        double yAvg = 1. / magnetsCount;
        sum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                sum += Math.sin(getAngleInRadians(lattice[i_row][i_col]));
            }
        }
        yAvg *= sum;
        return Math.sqrt(xAvg * xAvg + yAvg * yAvg);
    }

    private double calculateNearestNeighbourOrder(int magnetRow, int magnetCol) {
        int[][] lattice = latticeParametersImpl.lattice();
        ArrayList<Integer> neighboursStates = getNeighboursStates(lattice, magnetRow, magnetCol, 1);
        double onn = 1 / (magnetsCount * neighboursStates.size());
        double iSum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                double jSum = 0;
                for (int j = 0; j < neighboursStates.size(); j++) {
                    double alphaI = getAngleInRadians(lattice[i_row][i_col]);
                    double alphaJ = getAngleInRadians(neighboursStates.get(j));
                    jSum += Math.cos(alphaI - alphaJ);
                }
                iSum += jSum;
            }
        }
        onn *= iSum;
        return onn;
    }

    private double calculateP(ProbabilityFormula formula, double deltaE, double kBT) {
        if (formula == ProbabilityFormula.GLAUBER) {
            return Math.exp(-deltaE / kBT) / (1 + Math.exp(-deltaE / kBT));
        } else {
            if (deltaE > 0) {
                return Math.exp(-deltaE / kBT);
            } else {
                return 1;
            }
        }
    }

    private void changeMagnetState(int[][] lattice, int magnetRow, int magnetCol, int change) {
        int availableStates = latticeParametersImpl.states();
        int currentState = lattice[magnetRow][magnetCol];
        if (currentState + change >= availableStates || currentState + change < 0) {
            if (currentState + change >= availableStates) {
                lattice[magnetRow][magnetCol] = currentState + change - availableStates;
            } else {
                lattice[magnetRow][magnetCol] = availableStates + change;
            }
        } else {
            lattice[magnetRow][magnetCol] += change;
        }
    }

    public class LatticeParametersImpl implements LatticeParameters {

        private double _totalEnergy = 0;
        private double _orderParameter = 0;
        private double _nearestNeighbourOrder = 0;
        private int[][] _lattice;
        private int _states;

        @Override
        public double totalEnergy() {
            if (_totalEnergy == 0) {
                _totalEnergy = calculateTotalEnergy(_lattice);
            }
            return _totalEnergy;
        }

        @Override
        public double orderParameter() {
            if (_orderParameter == 0) {
                _orderParameter = calculateOrderParameter();
            }
            return _orderParameter;
        }

        @Override
        public double nearestNeighbourOrder() {
            if (_nearestNeighbourOrder == 0) {
                // Dlaczego ta metoda nearestNeighbourOrder nie przyjmuje argumentów określających, o który magnes konkretnie chodzi?
                _nearestNeighbourOrder = calculateNearestNeighbourOrder(0, 0);
            }
            return _nearestNeighbourOrder;
        }

        @Override
        public int[][] lattice() {
            return _lattice;
        }

        protected int states() {
            return _states;
        }

        protected int magnetState(int magnetRow, int magnetCol) {
            return _lattice[magnetRow][magnetCol];
        }

        protected void setLattice(int[][] lattice) {
            _lattice = new int[lattice.length][lattice.length];
            for (int i = 0; i < lattice.length; i++) {
                for (int j = 0; j < lattice.length; j++) {
                    _lattice[i][j] = lattice[i][j];
                }
            }
        }

        protected void setStates(int states) {
            _states = states;
        }

        protected void setTotalEnergy(double totalEnergy) {
            _totalEnergy = totalEnergy;
        }
    }
}