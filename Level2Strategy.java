import java.util.ArrayList;

public class Level2Strategy implements NeighboursStatesStrategy {

    @Override
    public ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol) {        
        ArrayList<Integer> neighboursStates = new ArrayList<Integer>();
        
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

        return neighboursStates;
    }
}
