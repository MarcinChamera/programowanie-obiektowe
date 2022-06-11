import java.util.ArrayList;

public class Level5Strategy implements NeighboursStatesStrategy {

    @Override
    public ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol) {        
        ArrayList<Integer> neighboursStates = new ArrayList<Integer>();
        
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

        return neighboursStates;
    }
}