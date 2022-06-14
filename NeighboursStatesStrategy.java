import java.util.ArrayList;

public interface NeighboursStatesStrategy {
    ArrayList<Integer> getNeighboursStates(int[][] lattice, int magnetRow, int magnetCol);
}
