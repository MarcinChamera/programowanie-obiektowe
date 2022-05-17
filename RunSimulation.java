import main.MCSimulation;
import main.Simulation.ProbabilityFormula;

import java.util.ArrayList;
import java.util.List;

/**
 * Symulacja kodu uruchamianego przez P. Oramusa.
 */
public class RunSimulation {
    public static void main(String[] args) {
        int steps = 5000;
        double TkB = 2.5;
        ProbabilityFormula formula = ProbabilityFormula.GLAUBER;
        int states = 2;
        int size = 10;
        double externaFieldAngle = 0;
        int[][] lattice = prepareArrangedLattice(size, 1);
        List<Double> parameters = new ArrayList<Double>();
        // Współczynnik oddziaływania z zewnętrznym polem
        Double Ce = 0.;
        parameters.add(Ce);
        // Współczynnik oddziaływania z sąsiadami, np. z najbliższymi to C1 (np. C1 = 1)
        //
        // Skąd wiadomo, ile współczynników trzeba dodać? 
        int Cn = 1;
        Double CValue = 1.;
        for(int i = 0; i < Cn; i++) {
            parameters.add(CValue);
        }

        MCSimulation mcSimulation = new MCSimulation();
        mcSimulation.setEnergyParameters(parameters, externaFieldAngle);
        mcSimulation.setLattice(lattice, states);
        mcSimulation.setProbabilityFormula(formula);
        mcSimulation.setTkB(TkB);

        // Czy program działa poprawnie
        double EStart = mcSimulation.getState().totalEnergy();
        System.out.println("Energia przed MC: " + Double.toString(EStart)); 
        mcSimulation.executeMCSteps(steps);
        double EEnd = mcSimulation.getState().totalEnergy();
        System.out.println("Energia po MC: " + Double.toString(EEnd)); 
        System.out.println("Energia na jeden magnes: " + Double.toString(EEnd / (size * size)));
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
}
