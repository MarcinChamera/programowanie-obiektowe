public class TotalEnergyDecorator implements TotalEnergy {
    protected TotalEnergy wrappee;

    public TotalEnergyDecorator(TotalEnergy totalEnergy) {
        wrappee = totalEnergy;
    }

    public double calculate(MCSimulation.LatticeParametersImpl latticeParametersImpl,  double currentResult) {
        return wrappee.calculate(latticeParametersImpl, currentResult);
    }
}
