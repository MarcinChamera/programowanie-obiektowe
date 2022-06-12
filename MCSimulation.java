import java.util.ArrayList;
import java.util.List;
import java.lang.Math;
import java.util.Random;

public class MCSimulation implements main.Simulation {

    private LatticeParametersImpl latticeParametersImpl = new LatticeParametersImpl();
    private ProbabilityFormula formula;
    private double TkB;
    private double Ce;
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
        latticeParametersImpl.setExternalFieldAngle(externaFieldAngle);
        latticeParametersImpl.setCn(parameters);
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
        double totalEnergy = latticeParametersImpl.totalEnergy();
        latticeParametersImpl.setTotalEnergy(totalEnergy);
        double acceptances = 0;
        double acceptanceRatio = 0;
        Random random; 
        for (int step = 0; step < steps; step++) {
            if (step > 0) {
                acceptanceRatio = acceptances / step;
            }
            double deltaE = 0;
            random = new Random();
            int magnetRowRandom = random.nextInt((int)Math.sqrt(magnetsCount));
            int magnetColRandom = random.nextInt((int)Math.sqrt(magnetsCount));
            int[][] newLattice = MCHelperSingleton.getInstance().generateLatticeCopy(latticeParametersImpl.lattice());
            if (acceptanceRatio > 0.5) {
                int magnetRowRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                int magnetColRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                while (magnetRowRandom == magnetRowRandom2 && magnetColRandom == magnetColRandom2) {
                    magnetRowRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                    magnetColRandom2 = random.nextInt((int)Math.sqrt(magnetsCount));
                }
                int magnetStateChange = random.nextBoolean() ? 1 : -1;
                int magnetStateChange2 = random.nextBoolean() ? 1 : -1;
                MCHelperSingleton.getInstance().changeMagnetState(newLattice, magnetRowRandom, magnetColRandom, magnetStateChange, latticeParametersImpl.states());
                MCHelperSingleton.getInstance().changeMagnetState(newLattice, magnetRowRandom2, magnetColRandom2, magnetStateChange2, latticeParametersImpl.states());
                deltaE = calculateTotalEnergy(newLattice) - totalEnergy;
            } else {
                int magnetStateChange = random.nextBoolean() ? 1 : -1;
                MCHelperSingleton.getInstance().changeMagnetState(newLattice, magnetRowRandom, magnetColRandom, magnetStateChange, latticeParametersImpl.states());
                deltaE = calculateEi(newLattice, magnetRowRandom, magnetColRandom) - calculateEi(latticeParametersImpl.lattice(), magnetRowRandom, magnetColRandom);
            }
            double R = random.nextDouble();
            double P = calculateP(formula, deltaE, TkB);
            if (R < P) {
                totalEnergy += deltaE;
                latticeParametersImpl.setTotalEnergy(totalEnergy);
                latticeParametersImpl.setLattice(newLattice);
                acceptances++;
            }
        }
    }

    private double calculateEi(int[][] lattice, int i_row, int i_col) {
        double Ei = 0;
        for (int n = 1; n < latticeParametersImpl.Cn().size(); n++) {
            if (latticeParametersImpl.Cn().get(n) == 0) {
                continue;
            }
            ArrayList<Integer> neighboursStates = MCHelperSingleton.getInstance().getNeighboursStates(lattice, i_row, i_col, n);
            for (int j = 0; j < neighboursStates.size(); j++) { 
                double alphaI = MCHelperSingleton.getInstance().getAngleInRadians(lattice[i_row][i_col], latticeParametersImpl.states());
                double alphaJ = MCHelperSingleton.getInstance().getAngleInRadians(neighboursStates.get(j), latticeParametersImpl.states());
                Ei -= latticeParametersImpl.Cn().get(n) * Math.cos(alphaI - alphaJ);
            }
        }
        return Ei;
    }

    private double calculateTotalEnergy(int[][] lattice) {
        TotalEnergy totalEnergy = new NoSubtractTotalEnergy();
        if (Ce != 0) {
            totalEnergy = new SubtractDecorator(totalEnergy);
        }
        return totalEnergy.calculate(latticeParametersImpl, -0.5);
    }

    private double calculateOrderParameter() {
        int[][] lattice = latticeParametersImpl.lattice();
        double xAvg = 1. / magnetsCount;
        double sum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                sum += Math.cos(MCHelperSingleton.getInstance().getAngleInRadians(lattice[i_row][i_col], latticeParametersImpl.states()));
            }
        }
        xAvg *= sum;
        double yAvg = 1. / magnetsCount;
        sum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                sum += Math.sin(MCHelperSingleton.getInstance().getAngleInRadians(lattice[i_row][i_col], latticeParametersImpl.states()));
            }
        }
        yAvg *= sum;
        return Math.sqrt(xAvg * xAvg + yAvg * yAvg);
    }

    private double calculateNearestNeighbourOrder() {
        int[][] lattice = latticeParametersImpl.lattice();
        ArrayList<Integer> neighboursStates = MCHelperSingleton.getInstance().getNeighboursStates(lattice, 0, 0, 1);
        double onn = 1. / (double)(magnetsCount * neighboursStates.size());
        double iSum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                double jSum = 0;
                neighboursStates = MCHelperSingleton.getInstance().getNeighboursStates(lattice, i_row, i_col, 1);
                for (int j = 0; j < neighboursStates.size(); j++) {
                    double alphaI = MCHelperSingleton.getInstance().getAngleInRadians(lattice[i_row][i_col], latticeParametersImpl.states());
                    double alphaJ = MCHelperSingleton.getInstance().getAngleInRadians(neighboursStates.get(j), latticeParametersImpl.states());
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

    public class LatticeParametersImpl implements LatticeParameters {

        private double _totalEnergy = 0;
        private int[][] _lattice;
        private int _states;
        private double _externalFieldAngle;
        private List<Double> _Cn;

        @Override
        public double totalEnergy() {
            if (_totalEnergy == 0) {
                _totalEnergy = calculateTotalEnergy(_lattice);
            }
            return _totalEnergy;
        }

        @Override
        public double orderParameter() {
            return calculateOrderParameter();
        }

        @Override
        public double nearestNeighbourOrder() {
            return calculateNearestNeighbourOrder();
        }

        @Override
        public int[][] lattice() {
            return _lattice;
        }

        protected int states() {
            return _states;
        }

        public double externalFieldAngle() {
            return _externalFieldAngle;
        }

        public List<Double> Cn() {
            return _Cn;
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

        protected void setExternalFieldAngle(double externalFieldAngle) {
            this._externalFieldAngle = externalFieldAngle;
        }

        protected void setCn(List<Double> Cn) {
            this._Cn = Cn;
        }
    }
}