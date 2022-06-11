import java.util.ArrayList;

public class Level1Strategy implements NeighboursStatesStrategy {

    @Override
    public ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol) {        
        ArrayList<Integer> neighboursStates = new ArrayList<Integer>();

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

        return neighboursStates;
    }   
}
