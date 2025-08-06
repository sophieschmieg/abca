package util.graph;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.vectors.Vector;
import util.MiscAlgorithms;
import util.graph.Graph.GraphBuilder;

class SudokoTest {

	private Graph<Vector<IntE>, IntE> sudokuGraph() {
		Integers z = Integers.z();
		GraphBuilder<Vector<IntE>, IntE> builder = new GraphBuilder<>(false);
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				builder.addVertex(new Vector<>(z.getInteger(i), z.getInteger(j)));
				for (int k = 0; k < 9; k++) {
					if (i != k) {
						builder.addEdge(new Vector<>(z.getInteger(i), z.getInteger(j)),
								new Vector<>(z.getInteger(k), z.getInteger(j)), z.zero());
					}
					if (j != k) {
						builder.addEdge(new Vector<>(z.getInteger(i), z.getInteger(j)),
								new Vector<>(z.getInteger(i), z.getInteger(k)), z.zero());
					}
					int x = 3 * MiscAlgorithms.DivRoundDown(i, 3) + (k % 3);
					int y = 3 * MiscAlgorithms.DivRoundDown(j, 3) + MiscAlgorithms.DivRoundDown(k, 3);
					if (i != x && j != y) {
						builder.addEdge(new Vector<>(z.getInteger(i), z.getInteger(j)),
								new Vector<>(z.getInteger(x), z.getInteger(y)), z.zero());
					}
				}
			}
		}
		return builder.build();
	}

	@Test
	void test() {
		Integers z = Integers.z();
		Graph<Vector<IntE>, IntE> graph = sudokuGraph();
		Set<IntE> colors = new TreeSet<>();
		for (int i = 0; i < 9; i++) {
			colors.add(z.getInteger(i + 1));
		}
		Map<Vector<IntE>, IntE> coloring = new TreeMap<>();
		coloring.put(new Vector<>(z.getInteger(0), z.getInteger(0)), z.getInteger(4));
		coloring.put(new Vector<>(z.getInteger(0), z.getInteger(1)), z.getInteger(1));
		coloring.put(new Vector<>(z.getInteger(0), z.getInteger(7)), z.getInteger(5));
		coloring.put(new Vector<>(z.getInteger(0), z.getInteger(8)), z.getInteger(8));
		coloring.put(new Vector<>(z.getInteger(1), z.getInteger(1)), z.getInteger(7));
		coloring.put(new Vector<>(z.getInteger(1), z.getInteger(8)), z.getInteger(6));
		coloring.put(new Vector<>(z.getInteger(2), z.getInteger(2)), z.getInteger(8));
		coloring.put(new Vector<>(z.getInteger(2), z.getInteger(4)), z.getInteger(1));
		coloring.put(new Vector<>(z.getInteger(2), z.getInteger(6)), z.getInteger(7));
		coloring.put(new Vector<>(z.getInteger(3), z.getInteger(5)), z.getInteger(4));
		coloring.put(new Vector<>(z.getInteger(3), z.getInteger(6)), z.getInteger(5));
		coloring.put(new Vector<>(z.getInteger(3), z.getInteger(7)), z.getInteger(9));
		coloring.put(new Vector<>(z.getInteger(3), z.getInteger(8)), z.getInteger(2));
		coloring.put(new Vector<>(z.getInteger(4), z.getInteger(3)), z.getInteger(8));
		coloring.put(new Vector<>(z.getInteger(4), z.getInteger(5)), z.getInteger(6));
		coloring.put(new Vector<>(z.getInteger(5), z.getInteger(0)), z.getInteger(9));
		coloring.put(new Vector<>(z.getInteger(5), z.getInteger(1)), z.getInteger(4));
		coloring.put(new Vector<>(z.getInteger(5), z.getInteger(2)), z.getInteger(7));
		coloring.put(new Vector<>(z.getInteger(5), z.getInteger(3)), z.getInteger(3));
		coloring.put(new Vector<>(z.getInteger(6), z.getInteger(2)), z.getInteger(5));
		coloring.put(new Vector<>(z.getInteger(6), z.getInteger(4)), z.getInteger(6));
		coloring.put(new Vector<>(z.getInteger(6), z.getInteger(6)), z.getInteger(3));
		coloring.put(new Vector<>(z.getInteger(7), z.getInteger(0)), z.getInteger(7));
		coloring.put(new Vector<>(z.getInteger(7), z.getInteger(7)), z.getInteger(6));
		coloring.put(new Vector<>(z.getInteger(8), z.getInteger(0)), z.getInteger(6));
		coloring.put(new Vector<>(z.getInteger(8), z.getInteger(1)), z.getInteger(3));
		coloring.put(new Vector<>(z.getInteger(8), z.getInteger(7)), z.getInteger(4));
		coloring.put(new Vector<>(z.getInteger(8), z.getInteger(8)), z.getInteger(5));
		Optional<Map<Vector<IntE>, IntE>> solved = graph.findVertexColoring(coloring, colors);
		assertTrue(solved.isPresent());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				System.out.print(solved.get().get(new Vector<>(z.getInteger(i), z.getInteger(j))));
				System.out.print(" ");
			}
			System.out.println();
		}
	}

}
