// import mathExt.Average;
// import mathExt.Energy;
// import lattice.LatticeState;
// import org.junit.jupiter.api.BeforeEach;
// import org.junit.jupiter.api.DisplayName;
// import org.junit.jupiter.params.ParameterizedTest;
// import org.junit.jupiter.params.provider.Arguments;
// import org.junit.jupiter.params.provider.MethodSource;
// import org.mockito.Mockito;

// import java.lang.reflect.Field;
// import java.util.Arrays;
// import java.util.LinkedList;
// import java.util.List;
// import java.util.stream.Stream;

// import static java.lang.Math.abs;
// import static main.Simulation.ProbabilityFormula.METROPOLIS;
// import static mathExt.MathExtended.angleInRad;
// import static org.junit.jupiter.api.Assertions.*;
// import static org.mockito.Mockito.*;

// public class MCSimulationTest {
//     MCSimulation mcSimulation;

//     static double DELTA = 0.0000001;

//     public static int[][] latticeTemplate = {
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//             {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

//     @BeforeEach
//     void setUp() {
// //        mcSimulation = new MCSimulation();
//         resetSingleton();
//         mcSimulation = Mockito.mock(MCSimulation.class);
// //        Energy.instance = new Energy();

// //        Mockito.when(mcSimulation.getState()).thenCallRealMethod(); // getter
//         doCallRealMethod().when(mcSimulation).setTkB(2.5); // setter
//         doCallRealMethod().when(mcSimulation).setProbabilityFormula(METROPOLIS);
//         mcSimulation.setTkB(2.5);
//         mcSimulation.setProbabilityFormula(METROPOLIS);

//         List<Double> energy = new LinkedList<>(Arrays.asList(0.0, 1.0));
//         doCallRealMethod().when(mcSimulation).setEnergyParameters(energy, 0);
//         mcSimulation.setEnergyParameters(energy, 0);


//         doCallRealMethod().when(mcSimulation).setLattice(latticeTemplate, 2);
//         mcSimulation.setLattice(latticeTemplate, 2);
//     }

//     public static void resetSingleton() {
//         Field field = null;
//         try {
//             field = Energy.class.getDeclaredField("instance");
//             field.setAccessible(true);
//             field.set(null, new Energy());
//         } catch (NoSuchFieldException | IllegalAccessException e) {
//             e.printStackTrace();
//         }
//     }

//     @ParameterizedTest
//     @MethodSource("initialTestProvider")
//     @DisplayName("Simple initial test")
//     void testInitial(int[][] lattice, int states, double nno, double te, double op) {
//         doCallRealMethod().when(mcSimulation).setLattice(lattice, states);
//         Mockito.when(mcSimulation.getState()).thenCallRealMethod();
//         mcSimulation.setLattice(lattice, states);

//         assertArrayEquals(mcSimulation.getState().lattice(), lattice);
//         assertEquals(mcSimulation.getState().nearestNeighbourOrder(), nno);
//         assertEquals(mcSimulation.getState().totalEnergy(), te);
//         assertEquals(mcSimulation.getState().orderParameter(), op);
//     }

//     static Stream<Arguments> initialTestProvider() {
//         int[][] lattice = {
//                 {1, 1, 1},
//                 {1, 1, 1},
//                 {1, 1, 1}};

//         return Stream.of(
//                 Arguments.of(lattice, 2, 1., -18., 1.)
//         );
//     }

//     @ParameterizedTest
//     @MethodSource("nearestNeighbourOrderTestProvider")
//     void testNearestNeighbourOrder(int[][] lattice, int states, double expectedResult) {
//         doCallRealMethod().when(mcSimulation).setLattice(lattice, states);
//         Mockito.when(mcSimulation.getState()).thenCallRealMethod();
//         mcSimulation.setLattice(lattice, states);

//         double result = mcSimulation.getState().nearestNeighbourOrder();

//         assertTrue(result >= -1);
//         assertTrue(result <= 1);
//         assertTrue(abs(expectedResult - result) < DELTA);
//     }

//     static Stream<Arguments> nearestNeighbourOrderTestProvider() {

//         int[][] lattice = {
//                 {1, 1, 1},
//                 {1, 1, 1},
//                 {1, 1, 1}};

//         int[][] lattice2 = {
//                 {0, 1, 0},
//                 {1, 0, 1},
//                 {0, 1, 0}};

//         int[][] lattice3 = {
//                 {0, 1, 0, 1},
//                 {1, 0, 1, 0},
//                 {0, 1, 0, 1},
//                 {1, 0, 1, 0}};

//         int[][] lattice4 = {
//                 {0, 0, 0},
//                 {0, 0, 0},
//                 {0, 0, 0}};

//         return Stream.of(
//                 Arguments.of(latticeTemplate, 2, 1.),
//                 Arguments.of(latticeTemplate, 3, 1.),
//                 Arguments.of(lattice, 2, 1.),
//                 Arguments.of(lattice, 3, 1.),
//                 Arguments.of(lattice2, 2, -0.33333333333),
//                 Arguments.of(lattice3, 2, -1.),
//                 Arguments.of(lattice3, 3, -0.5),
//                 Arguments.of(lattice4, 2, 1.),
//                 Arguments.of(lattice4, 3, 1.)
//         );
//     }
// }