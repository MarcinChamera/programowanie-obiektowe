import main.Simulation.ProbabilityFormula;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class SimulationTests {
    public static void main(String[] args) {
        RandomLatticeTest();
        ArrangedLatticeTest();
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

    public static void RandomLatticeTest() {
        int steps = 1000;
        double TkB = 2.5;
        ProbabilityFormula formula = ProbabilityFormula.GLAUBER;
        int states = 2;
        int size = 32;
        double externaFieldAngle = 0;
        // int[][] lattice = prepareArrangedLattice(size, 1);
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
        System.out.println("Energia na jeden magnes: " + Double.toString(afterOneMagnesEnergy));
        double afterOrderParameter =  mcSimulation.getState().orderParameter();
        System.out.println("Order parameter: " + Double.toString(afterOrderParameter));
        double afterNNOrderParameter =  mcSimulation.getState().nearestNeighbourOrder();
        System.out.println("NNOrder parameter: " + Double.toString(afterNNOrderParameter));

        // double energyDifferenceThreshold = 10;
        // double oneMagnesTrueEnergy = 0;
        // double trueOrderParameter = 0;
        // double trueNNOrderParameter = 0;
        // double oneMagnesBeforeEnergyDifference = oneMagnesTrueEnergy - beforeOneMagnesEnergy;
        // double oneMagnesAfterEnergyDifference = oneMagnesTrueEnergy - afterOneMagnesEnergy;

        // System.out.println("\n");

        // if (Math.abs(EStart - EEnd) > energyDifferenceThreshold) System.out.println("Roznica energii przed i po wieksza o " + Double.toString(Math.abs(EStart - EEnd)));
        // if (oneMagnesTrueEnergy) System.out.println("Energia na jeden magnes przed symulacja rozni sie o wiecej niz " + Double.toString(oneMagnesBeforeEnergyDifference));
        // if (beforeOrderParameter != trueOrderParameter) System.out.println("Parametr uporzadkowanie przed symulacja rozni sie od poprawnego " + Double.toString(beforeOrderParameter));
        // if (beforeNNOrderParameter != trueNNOrderParameter) System.out.println("Parametr uporzadkowanie bliskiego zasiegu przed symulacja rozni sie od poprawnego " + Double.toString(beforeNNOrderParameter));
        // if (Math.abs(oneMagnesAfterEnergyDifference) > 0.1) System.out.println("Energia na jeden magnes po symulacji rozni sie o wiecej niz " + Double.toString(oneMagnesAfterEnergyDifference));
        // if (afterOrderParameter != trueOrderParameter) System.out.println("Parametr uporzadkowanie po symulacji rozni sie od poprawnego o " + Double.toString(Math.abs(afterOrderParameter - trueOrderParameter)));
        // if (afterNNOrderParameter != trueNNOrderParameter) System.out.println("Parametr uporzadkowanie bliskiego zasiegu po symulacji rozni sie od poprawnego o " + Double.toString(Math.abs(afterNNOrderParameter - trueNNOrderParameter)));
        System.out.println("=======================");
    }

    public static void ArrangedLatticeTest() {
        int steps = 1000;
        double TkB = 2.5;
        ProbabilityFormula formula = ProbabilityFormula.GLAUBER;
        int states = 2;
        int size = 32;
        double externaFieldAngle = 0;
        int[][] lattice = prepareArrangedLattice(size, 1);
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
        System.out.println("Energia na jeden magnes (Etot / (N*N)): " + Double.toString(afterOneMagnesEnergy));
        double afterOrderParameter =  mcSimulation.getState().orderParameter();
        System.out.println("Order parameter: " + Double.toString(afterOrderParameter));
        double afterNNOrderParameter =  mcSimulation.getState().nearestNeighbourOrder();
        System.out.println("NNOrder parameter: " + Double.toString(afterNNOrderParameter));

        // double energyDifferenceThreshold = 10;
        // double oneMagnesTrueEnergy = 2;
        // double trueOrderParameter = 1;
        // double trueNNOrderParameter = 1;
        // double oneMagnesBeforeEnergyDifference = oneMagnesTrueEnergy - beforeOneMagnesEnergy;
        // double oneMagnesAfterEnergyDifference = oneMagnesTrueEnergy - afterOneMagnesEnergy;

        // System.out.println("\n");

        // if (Math.abs(EStart - EEnd) > energyDifferenceThreshold) System.out.println("Roznica energii przed i po wieksza o " + Double.toString(Math.abs(EStart - EEnd)));
        // if (Math.abs(oneMagnesBeforeEnergyDifference) > 0.1) System.out.println("Energia na jeden magnes przed symulacja rozni sie o wiecej niz " + Double.toString(oneMagnesBeforeEnergyDifference));
        // if (beforeOrderParameter != trueOrderParameter) System.out.println("Parametr uporzadkowanie przed symulacja rozni sie od poprawnego " + Double.toString(beforeOrderParameter));
        // if (beforeNNOrderParameter != trueNNOrderParameter) System.out.println("Parametr uporzadkowanie bliskiego zasiegu przed symulacja rozni sie od poprawnego " + Double.toString(beforeNNOrderParameter));
        // if (Math.abs(oneMagnesAfterEnergyDifference) > 0.1) System.out.println("Energia na jeden magnes po symulacji rozni sie o wiecej niz " + Double.toString(oneMagnesAfterEnergyDifference));
        // if (afterOrderParameter != trueOrderParameter) System.out.println("Parametr uporzadkowanie po symulacji rozni sie od poprawnego o " + Double.toString(Math.abs(afterOrderParameter - trueOrderParameter)));
        // if (afterNNOrderParameter != trueNNOrderParameter) System.out.println("Parametr uporzadkowanie bliskiego zasiegu po symulacji rozni sie od poprawnego o " + Double.toString(Math.abs(afterNNOrderParameter - trueNNOrderParameter)));
        System.out.println("=======================");
    }
}
