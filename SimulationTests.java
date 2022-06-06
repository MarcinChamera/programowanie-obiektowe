import main.Simulation.ProbabilityFormula;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

class OnnParameterTest {
    double DELTA = 0.0000001;

    int[][] latticeTemplate = {
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

    int[][] lattice = {
        {1, 1, 1},
        {1, 1, 1},
        {1, 1, 1}};

    int[][] lattice2 = {
            {0, 1, 0},
            {1, 0, 1},
            {0, 1, 0}};

    int[][] lattice3 = {
            {0, 1, 0, 1},
            {1, 0, 1, 0},
            {0, 1, 0, 1},
            {1, 0, 1, 0}};

    int[][] lattice4 = {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}};

    public void run() {
        singleOnnTest(latticeTemplate, 2, 1., "Onn Test 1: ");
        singleOnnTest(latticeTemplate, 3, 1., "Onn Test 2: ");
        singleOnnTest(lattice, 2, 1., "Onn Test 3: ");
        singleOnnTest(lattice, 3, 1., "Onn Test 4: ");
        singleOnnTest(lattice2, 2, -0.33333333333, "Onn Test 5: ");
        singleOnnTest(lattice3, 2, -1., "Onn Test 6: ");
        singleOnnTest(lattice3, 3, -0.5, "Onn Test 7: ");
        singleOnnTest(lattice4, 2, 1., "Onn Test 8: ");
        singleOnnTest(lattice4, 3, 1., "Onn Test 9: ");
    }

    private void singleOnnTest(int[][] lattice, int states, double expectedResult, String message) {
        MCSimulation mcSimulation = new MCSimulation();
        mcSimulation.setLattice(lattice, states);
        double result = mcSimulation.getState().nearestNeighbourOrder();
        boolean passed = Math.abs(expectedResult - result) < DELTA;
        if (passed) {
            System.out.println(message + passed);
        } else {
            System.out.println(message + passed + " expected: " + expectedResult + " result: " + result);
        }
    }
}

public class SimulationTests {
    private static ArrayList<Double> endEnergyPerMagnes;
    static double DELTA = 0.0000001;

    public static void main(String[] args) {
        // runArrangedLatticeTestNTimes(10);
        // randomLatticeTest();
        new OnnParameterTest().run();
    }

    private static void runArrangedLatticeTestNTimes(int n) {
        endEnergyPerMagnes = new ArrayList<>();
        double sum = 0;
        for (int i = 0; i < n; i++) {
            arrangedLatticeTest();
            sum += endEnergyPerMagnes.get(endEnergyPerMagnes.size()-1);
        }
        System.out.println("Srednia koncowa energia jednego magnesu = " + Double.toString(sum / endEnergyPerMagnes.size()));
        System.out.println("=======================");
    }

    private static int[][] prepareRandomLattice(int size, int states) {
        int[][] lattice = new int[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int rand = new Random().nextInt(states);
                lattice[i][j] = rand;
            }
        }
        return lattice;
    }

    private static int[][] prepareArrangedLattice(int size, int state) {
        int[][] lattice = new int[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                lattice[i][j] = state;
            }
        }
        return lattice;
    }

    // tylko dla dwoch stanow na razie
    // tylko dla parzystych wymiarow na razie
    private static int[][] prepareAlternatingLattice(int size, int states) {
        int[][] lattice = new int[size][size];
        int currentState = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                lattice[i][j] = currentState;
                currentState++;
                currentState = currentState % 2;
            }
        }
        return lattice;
    }

    public static void randomLatticeTest() {
        int steps = 10000;
        double TkB = 2.5;
        ProbabilityFormula formula = ProbabilityFormula.METROPOLIS;
        int states = 2;
        int size = 32;
        double externaFieldAngle = 0;
        int[][] lattice = prepareRandomLattice(size, states);
        List<Double> parameters = new ArrayList<Double>();
        // Współczynnik oddziaływania z zewnętrznym polem
        Double Ce = 0.;
        parameters.add(Ce);
        // Współczynnik oddziaływania z sąsiadami, np. z najbliższymi to C1 (np. C1 = 1)
        Double C1 = 1.;
        parameters.add(C1);

        MCSimulation mcSimulation = new MCSimulation();
        mcSimulation.setEnergyParameters(parameters, externaFieldAngle);
        mcSimulation.setLattice(lattice, states);
        mcSimulation.setProbabilityFormula(formula);
        mcSimulation.setTkB(TkB);

        System.out.println("Test tablicy losowym uporzadkowaniem\n");
        double EStart = mcSimulation.getState().totalEnergy();
        System.out.println("Energia przed MC: " + Double.toString(EStart));
        double beforeOneMagnesEnergy =  EStart / (size * size);
        System.out.println("Energia na jeden magnes: " + Double.toString(beforeOneMagnesEnergy));
        double beforeOrderParameter =  mcSimulation.getState().orderParameter();
        System.out.println("Order parameter: " + Double.toString(beforeOrderParameter));
        double beforeNNOrderParameter =  mcSimulation.getState().nearestNeighbourOrder();
        System.out.println("NNOrder parameter: " + Double.toString(beforeNNOrderParameter));

        mcSimulation.executeMCSteps(steps);

        double EEnd = mcSimulation.getState().totalEnergy();
        System.out.println("\nEnergia po MC: " + Double.toString(EEnd)); 
        double afterOneMagnesEnergy =  EEnd / (size * size);
        endEnergyPerMagnes.add(afterOneMagnesEnergy);
        System.out.println("Energia na jeden magnes: " + Double.toString(afterOneMagnesEnergy));
        double afterOrderParameter =  mcSimulation.getState().orderParameter();
        System.out.println("Order parameter: " + Double.toString(afterOrderParameter));
        double afterNNOrderParameter =  mcSimulation.getState().nearestNeighbourOrder();
        System.out.println("NNOrder parameter: " + Double.toString(afterNNOrderParameter));
        System.out.println("=======================");
    }

    public static void arrangedLatticeTest() {
        int steps = 10000;
        double TkB = 2.5;
        ProbabilityFormula formula = ProbabilityFormula.METROPOLIS;
        int states = 2;
        int size = 32;
        double externaFieldAngle = 0;
        int[][] lattice = prepareArrangedLattice(size, 0);
        List<Double> parameters = new ArrayList<Double>();
        // Współczynnik oddziaływania z zewnętrznym polem
        Double Ce = 0.;
        parameters.add(Ce);
        // Współczynnik oddziaływania z sąsiadami, np. z najbliższymi to C1 (np. C1 = 1)
        Double C1 = 1.;
        parameters.add(C1);

        MCSimulation mcSimulation = new MCSimulation();
        mcSimulation.setEnergyParameters(parameters, externaFieldAngle);
        mcSimulation.setLattice(lattice, states);
        mcSimulation.setProbabilityFormula(formula);
        mcSimulation.setTkB(TkB);

        System.out.println("Test uporzadkowanej tablicy\n");
        double EStart = mcSimulation.getState().totalEnergy();
        System.out.println("Energia przed MC: " + Double.toString(EStart));
        double beforeOneMagnesEnergy =  EStart / (size * size);
        System.out.println("Energia na jeden magnes (Etot / (N*N)): " + Double.toString(beforeOneMagnesEnergy));
        double beforeOrderParameter =  mcSimulation.getState().orderParameter();
        System.out.println("Order parameter: " + Double.toString(beforeOrderParameter));
        double beforeNNOrderParameter =  mcSimulation.getState().nearestNeighbourOrder();
        System.out.println("NNOrder parameter: " + Double.toString(beforeNNOrderParameter));

        mcSimulation.executeMCSteps(steps);

        double EEnd = mcSimulation.getState().totalEnergy();
        System.out.println("\nEnergia po MC: " + Double.toString(EEnd)); 
        double afterOneMagnesEnergy =  EEnd / (size * size);
        endEnergyPerMagnes.add(afterOneMagnesEnergy);
        System.out.println("Energia na jeden magnes (Etot / (N*N)): " + Double.toString(afterOneMagnesEnergy));
        double afterOrderParameter =  mcSimulation.getState().orderParameter();
        System.out.println("Order parameter: " + Double.toString(afterOrderParameter));
        double afterNNOrderParameter =  mcSimulation.getState().nearestNeighbourOrder();
        System.out.println("NNOrder parameter: " + Double.toString(afterNNOrderParameter));
        System.out.println("=======================");
    }
}
