public class TotalEnergyDecorator implements TotalEnergy {
    private TotalEnergy wrappee;

    public TotalEnergyDecorator(TotalEnergy totalEnergy) {
        wrappee = totalEnergy;
    }

    @Override
    public double calculate(MCSimulation.LatticeParametersImpl latticeParametersImpl,  double currentResult) {
        return wrappee.calculate(latticeParametersImpl, currentResult);
    }
}
