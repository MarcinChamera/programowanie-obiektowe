import java.util.ArrayList;

public class Level3Strategy implements NeighboursStatesStrategy {

    @Override
    public ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol) {        
        ArrayList<Integer> neighboursStates = new ArrayList<Integer>();
        
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

        return neighboursStates;
    }
}
