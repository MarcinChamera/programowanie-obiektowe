import java.util.ArrayList;

public class Level4Strategy implements NeighboursStatesStrategy {

    @Override
    public ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol) {        
        ArrayList<Integer> neighboursStates = new ArrayList<Integer>();
        
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

        return neighboursStates;
    }
}
