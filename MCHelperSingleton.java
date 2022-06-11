import java.util.ArrayList;

public class MCHelperSingleton {
    private NeighboursStatesStrategy neighboursStatesStrategy;
    private static MCHelperSingleton instance;
    private MCHelperSingleton() {};
    public static MCHelperSingleton getInstance() {
        if (instance == null) {
            instance = new MCHelperSingleton();
        }
        return instance;
    }

    public ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol, int level) {
        if (level == 1) {
            neighboursStatesStrategy = new Level1Strategy();
        } else if (level == 2) {
            neighboursStatesStrategy = new Level2Strategy();
        } else if (level == 3) {
            neighboursStatesStrategy = new Level3Strategy();
        } else if (level == 4) {
            neighboursStatesStrategy = new Level4Strategy();
        } else if (level == 5) {
            neighboursStatesStrategy = new Level5Strategy();
        }
        return neighboursStatesStrategy.getNeighboursStates(lattice, magnetRow, magnetCol);
    }

    public double getAngleInRadians(int magnetState, int states) {
        return 2 * Math.PI * magnetState / states;
    }

    public int[][] generateLatticeCopy(int[][] lattice) {
        int[][] latticeCopy = new int[lattice.length][lattice.length];
        for (int i = 0; i < lattice.length; i++) {
            for (int j = 0; j < lattice.length; j++) {
                latticeCopy[i][j] = lattice[i][j];
            }
        }
        return latticeCopy;
    }

    public void changeMagnetState(int[][] lattice, int magnetRow, int magnetCol, int change, int states) {
        int availableStates = states;
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
}
