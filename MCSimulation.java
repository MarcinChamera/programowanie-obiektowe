import java.util.ArrayList;
import java.util.List;
import java.lang.Math;
import java.util.Random;

public class MCSimulation implements main.Simulation {

    private LatticeParametersImpl latticeParametersImpl = new LatticeParametersImpl();
    private ProbabilityFormula formula;
    private double TkB;
    private double Ce;
    private List<Double> Cn;
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
        Cn = parameters;
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
            int[][] newLattice = generateLatticeCopy(latticeParametersImpl.lattice());
            if (acceptanceRatio > 0.5) {
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
                deltaE = calculateTotalEnergy(newLattice) - totalEnergy;
            } else {
                int magnetStateChange = random.nextBoolean() ? 1 : -1;
                changeMagnetState(newLattice, magnetRowRandom, magnetColRandom, magnetStateChange);
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
        for (int n = 1; n < Cn.size(); n++) {
            ArrayList<Integer> neighboursStates = getNeighboursStates(lattice, i_row, i_col, n);
            for (int j = 0; j < neighboursStates.size(); j++) { 
                double alphaI = getAngleInRadians(lattice[i_row][i_col]);
                double alphaJ = getAngleInRadians(neighboursStates.get(j));
                Ei -= Cn.get(n) * Math.cos(alphaI - alphaJ);
            }
        }
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
                        double alphaI = getAngleInRadians(lattice[i_row][i_col]);
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
                double alphaI = getAngleInRadians(lattice[i_row][i_col]);
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
            // gora
            if (magnetRow > 0) neighboursStates.add(lattice[magnetRow - 1][magnetCol]);
            else neighboursStates.add(lattice[lattice.length - 1][magnetCol]);

            // dol
            if (magnetRow < lattice.length - 1) neighboursStates.add(lattice[magnetRow + 1][magnetCol]);
            else neighboursStates.add(lattice[0][magnetCol]);

            // lewo
            if (magnetCol > 0) neighboursStates.add(lattice[magnetRow][magnetCol - 1]);
            else neighboursStates.add(lattice[magnetRow][lattice.length - 1]);

            // prawo
            if (magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow][magnetCol + 1]);
            else neighboursStates.add(lattice[magnetRow][0]);

        } else if (level == 2) {
            // gora-lewo
            if (magnetRow > 0 && magnetCol > 0) neighboursStates.add(lattice[magnetRow - 1][magnetCol - 1]);
            else if (magnetRow == 0 && magnetCol > 0) neighboursStates.add(lattice[lattice.length - 1][magnetCol - 1]);
            else if (magnetRow > 0 && magnetCol == 0) neighboursStates.add(lattice[magnetRow - 1][lattice.length - 1]);
            else if (magnetRow == 0 && magnetCol == 0) neighboursStates.add(lattice[lattice.length - 1][lattice.length - 1]);
            
            // dol-lewo
            if (magnetRow < lattice.length - 1 && magnetCol > 0) neighboursStates.add(lattice[magnetRow + 1][magnetCol - 1]);
            else if (magnetRow < lattice.length - 1 && magnetCol == 0) neighboursStates.add(lattice[magnetRow + 1][lattice.length - 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol == 0) neighboursStates.add(lattice[0][lattice.length - 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol > 0) neighboursStates.add(lattice[0][magnetCol - 1]);

            // gora-prawo
            if (magnetRow > 0 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow - 1][magnetCol + 1]);
            else if (magnetRow > 0 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow - 1][0]);
            else if (magnetRow == 0 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[lattice.length - 1][0]);
            else if (magnetRow == 0 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[lattice.length - 1][magnetCol + 1]);

            // dol-prawo
            if (magnetRow < lattice.length - 1 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow + 1][magnetCol + 1]);
            else if (magnetRow < lattice.length - 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow + 1][0]);
            else if (magnetRow == lattice.length - 1 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[0][magnetCol + 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[0][0]);

        } else if (level == 3) {
            // gora
            if (magnetRow > 1) neighboursStates.add(lattice[magnetRow - 2][magnetCol]);
            else if (magnetRow == 1) neighboursStates.add(lattice[lattice.length - 1][magnetCol]);
            else neighboursStates.add(lattice[lattice.length - 2][magnetCol]);

            // dol
            if (magnetRow < lattice.length - 2) neighboursStates.add(lattice[magnetRow + 2][magnetCol]);
            else if (magnetRow == lattice.length - 2) neighboursStates.add(lattice[0][magnetCol]);
            else neighboursStates.add(lattice[1][magnetCol]);

            // lewo
            if (magnetCol > 1) neighboursStates.add(lattice[magnetRow][magnetCol - 2]);
            else if (magnetCol == 1) neighboursStates.add(lattice[magnetRow][lattice.length - 1]);
            else neighboursStates.add(lattice[magnetRow][lattice.length - 2]);
            
            // prawo
            if (magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow][magnetCol + 2]);
            else if (magnetCol == lattice.length - 2) neighboursStates.add(lattice[magnetRow][0]);
            else neighboursStates.add(lattice[magnetRow][1]);

        } else if (level == 4) {
            // 1 do gory, 2 w lewo
            if (magnetRow > 0 && magnetCol > 1) neighboursStates.add(lattice[magnetRow - 1][magnetCol - 2]);
            else if (magnetRow > 0 && magnetCol == 1) neighboursStates.add(lattice[magnetRow - 1][lattice.length - 1]);
            else if (magnetRow > 0 && magnetCol == 0) neighboursStates.add(lattice[magnetRow - 1][lattice.length - 2]);
            else if (magnetRow == 0 && magnetCol > 1) neighboursStates.add(lattice[lattice.length - 1][magnetCol - 2]);
            else if (magnetRow == 0 && magnetCol == 1) neighboursStates.add(lattice[lattice.length - 1][lattice.length - 1]);
            else if (magnetRow == 0 && magnetCol == 0) neighboursStates.add(lattice[lattice.length - 1][lattice.length - 2]);

            // 1 w dol, 2 w lewo
            if (magnetRow < lattice.length - 1 && magnetCol > 1) neighboursStates.add(lattice[magnetRow + 1][magnetCol - 2]);
            else if (magnetRow < lattice.length - 1 && magnetCol == 1) neighboursStates.add(lattice[magnetRow + 1][lattice.length - 1]);
            else if (magnetRow < lattice.length - 1 && magnetCol == 0) neighboursStates.add(lattice[magnetRow + 1][lattice.length - 2]);
            else if (magnetRow == lattice.length - 1 && magnetCol > 1) neighboursStates.add(lattice[0][magnetCol - 2]);
            else if (magnetRow == lattice.length - 1 && magnetCol == 1) neighboursStates.add(lattice[0][lattice.length - 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol == 0) neighboursStates.add(lattice[0][lattice.length - 2]);

            // 1 w gore, 2 w prawo
            if (magnetRow > 0 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow - 1][magnetCol + 2]);
            else if (magnetRow > 0 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[magnetRow - 1][0]);
            else if (magnetRow > 0 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow - 1][1]);
            else if (magnetRow == 0 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[lattice.length - 1][magnetCol + 2]);
            else if (magnetRow == 0 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[lattice.length - 1][0]);
            else if (magnetRow == 0 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[lattice.length - 1][1]);

            // 1 w dol, 2 w prawo
            if (magnetRow < lattice.length - 1 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow + 1][magnetCol + 2]);
            else if (magnetRow < lattice.length - 1 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[magnetRow + 1][0]);
            else if (magnetRow < lattice.length - 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow + 1][1]);
            else if (magnetRow == lattice.length - 1 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[0][magnetCol + 2]);
            else if (magnetRow == lattice.length - 1 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[0][0]);
            else if (magnetRow == lattice.length - 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[0][1]);

            // 2 w gore, 1 w lewo
            if (magnetRow > 1 && magnetCol > 0) neighboursStates.add(lattice[magnetRow - 2][magnetCol - 1]);
            else if (magnetRow == 1 && magnetCol > 0) neighboursStates.add(lattice[lattice.length - 1][magnetCol - 1]);
            else if (magnetRow == 0 && magnetCol > 0) neighboursStates.add(lattice[lattice.length - 2][magnetCol - 1]);
            else if (magnetRow > 1 && magnetCol == 0) neighboursStates.add(lattice[magnetRow - 2][lattice.length - 1]);
            else if (magnetRow == 1 && magnetCol == 0) neighboursStates.add(lattice[lattice.length - 1][lattice.length - 1]);
            else if (magnetRow == 0 && magnetCol == 0) neighboursStates.add(lattice[lattice.length - 2][lattice.length - 1]);

            // 2 w dol, 1 w lewo
            if (magnetRow < lattice.length - 2 && magnetCol > 0) neighboursStates.add(lattice[magnetRow + 2][magnetCol - 1]);
            else if (magnetRow == lattice.length - 2 && magnetCol > 0) neighboursStates.add(lattice[0][magnetCol - 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol > 0) neighboursStates.add(lattice[1][magnetCol - 1]);
            else if (magnetRow < lattice.length - 2 && magnetCol == 0) neighboursStates.add(lattice[magnetRow + 2][lattice.length - 1]);
            else if (magnetRow == lattice.length - 2 && magnetCol == 0) neighboursStates.add(lattice[0][lattice.length - 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol == 0) neighboursStates.add(lattice[1][lattice.length - 1]);

            // 2 w gore, 1 w prawo
            if (magnetRow > 1 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow - 2][magnetCol + 1]);
            else if (magnetRow == 1 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[lattice.length - 1][magnetCol + 1]);
            else if (magnetRow == 0 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[lattice.length - 2][magnetCol + 1]);
            else if (magnetRow > 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow - 2][0]);
            else if (magnetRow == 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[lattice.length - 1][0]);
            else if (magnetRow == 0 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[lattice.length - 2][0]);

            // 2 w dol, 1 w prawo
            if (magnetRow < lattice.length - 2 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[magnetRow + 2][magnetCol + 1]);
            else if (magnetRow == lattice.length - 2 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[0][magnetCol + 1]);
            else if (magnetRow == lattice.length - 1 && magnetCol < lattice.length - 1) neighboursStates.add(lattice[1][magnetCol + 1]);
            else if (magnetRow < lattice.length - 2 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow + 2][0]);
            else if (magnetRow == lattice.length - 2 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[0][0]);
            else if (magnetRow == lattice.length - 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[1][0]);

        } else if (level == 5) {
            // gora-lewo
            if (magnetRow > 1 && magnetCol > 1) neighboursStates.add(lattice[magnetRow - 2][magnetCol - 2]);
            else if (magnetRow == 1 && magnetCol > 1) neighboursStates.add(lattice[lattice.length - 1][magnetCol - 2]);
            else if (magnetRow == 0 && magnetCol > 1) neighboursStates.add(lattice[lattice.length - 2][magnetCol - 2]);
            else if (magnetRow > 1 && magnetCol == 1) neighboursStates.add(lattice[magnetRow - 2][lattice.length - 1]);
            else if (magnetRow == 1 && magnetCol == 1) neighboursStates.add(lattice[lattice.length - 1][lattice.length - 1]);
            else if (magnetRow == 0 && magnetCol == 1) neighboursStates.add(lattice[lattice.length - 2][lattice.length - 1]);
            else if (magnetRow > 1 && magnetCol == 0) neighboursStates.add(lattice[magnetRow - 2][lattice.length - 2]);
            else if (magnetRow == 1 && magnetCol == 0) neighboursStates.add(lattice[lattice.length - 1][lattice.length - 2]);
            else if (magnetRow == 0 && magnetCol == 0) neighboursStates.add(lattice[lattice.length - 2][lattice.length - 2]);
            
            // dol-lewo
            if (magnetRow < lattice.length - 2 && magnetCol > 1) neighboursStates.add(lattice[magnetRow + 2][magnetCol - 2]);
            if (magnetRow == lattice.length - 2 && magnetCol > 1) neighboursStates.add(lattice[0][magnetCol - 2]);
            if (magnetRow == lattice.length - 1 && magnetCol > 1) neighboursStates.add(lattice[1][magnetCol - 2]);
            if (magnetRow < lattice.length - 2 && magnetCol == 1) neighboursStates.add(lattice[magnetRow + 2][lattice.length - 1]);
            if (magnetRow == lattice.length - 2 && magnetCol == 1) neighboursStates.add(lattice[0][lattice.length - 1]);
            if (magnetRow == lattice.length - 1 && magnetCol == 1) neighboursStates.add(lattice[1][lattice.length - 1]);
            if (magnetRow < lattice.length - 2 && magnetCol == 0) neighboursStates.add(lattice[magnetRow + 2][lattice.length - 2]);
            if (magnetRow == lattice.length - 2 && magnetCol == 0) neighboursStates.add(lattice[0][lattice.length - 2]);
            if (magnetRow == lattice.length - 1 && magnetCol == 0) neighboursStates.add(lattice[1][lattice.length - 2]);

            // gora-prawo
            if (magnetRow > 1 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow - 2][magnetCol + 2]);
            if (magnetRow == 1 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[lattice.length - 1][magnetCol + 2]);
            if (magnetRow == 0 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[lattice.length - 2][magnetCol + 2]);
            if (magnetRow > 1 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[magnetRow - 2][0]);
            if (magnetRow == 1 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[lattice.length - 1][0]);
            if (magnetRow == 0 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[lattice.length - 2][0]);
            if (magnetRow > 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow - 2][1]);
            if (magnetRow == 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[lattice.length - 1][1]);
            if (magnetRow == 0 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[lattice.length - 2][1]);

            // dol-prawo
            if (magnetRow < lattice.length - 2 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[magnetRow + 2][magnetCol + 2]);
            if (magnetRow == lattice.length - 2 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[0][magnetCol + 2]);
            if (magnetRow == lattice.length - 1 && magnetCol < lattice.length - 2) neighboursStates.add(lattice[1][magnetCol + 2]);
            if (magnetRow < lattice.length - 2 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[magnetRow + 2][0]);
            if (magnetRow == lattice.length - 2 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[0][0]);
            if (magnetRow == lattice.length - 1 && magnetCol == lattice.length - 2) neighboursStates.add(lattice[1][0]);
            if (magnetRow < lattice.length - 2 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[magnetRow + 2][1]);
            if (magnetRow == lattice.length - 2 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[0][1]);
            if (magnetRow == lattice.length - 1 && magnetCol == lattice.length - 1) neighboursStates.add(lattice[1][1]);
        }
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

    private double calculateNearestNeighbourOrder() {
        int[][] lattice = latticeParametersImpl.lattice();
        ArrayList<Integer> neighboursStates = getNeighboursStates(lattice, 0, 0, 1);
        double onn = 1. / (double)(magnetsCount * neighboursStates.size());
        double iSum = 0;
        for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
            for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                double jSum = 0;
                neighboursStates = getNeighboursStates(lattice, i_row, i_col, 1);
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