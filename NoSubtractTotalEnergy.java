import java.util.ArrayList;

public class NoSubtractTotalEnergy implements TotalEnergy {

    public double calculate(MCSimulation.LatticeParametersImpl latticeParametersImpl,  double currentResult) {
        int[][] lattice = latticeParametersImpl.lattice();
        int magnetsCount = lattice.length * lattice.length;
        double nSum = 0;
        for (int n = 1; n < latticeParametersImpl.Cn().size(); n++) {
            if (latticeParametersImpl.Cn().get(n) == 0) {
                continue;
            }
            double iSum = 0;
            for (int i_row = 0; i_row < Math.sqrt(magnetsCount); i_row++) {
                for (int i_col = 0; i_col < Math.sqrt(magnetsCount); i_col++) {
                    ArrayList<Integer> neighboursStates = MCHelperSingleton.getInstance().getNeighboursStates(lattice, i_row, i_col, n);
                    double jSum = 0;
                    for (int j = 0; j < neighboursStates.size(); j++) {
                        double alphaI = MCHelperSingleton.getInstance().getAngleInRadians(lattice[i_row][i_col], latticeParametersImpl.states());
                        double alphaJ = MCHelperSingleton.getInstance().getAngleInRadians(neighboursStates.get(j), latticeParametersImpl.states());
                        jSum +=  Math.cos(alphaI - alphaJ);
                    }
                    iSum += jSum;
                }
            }
            iSum *= latticeParametersImpl.Cn().get(n);
            nSum += iSum;
        }
        currentResult *= nSum;
        return currentResult;
    }
}
