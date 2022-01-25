package util.graph;

import java.util.List;

import org.junit.jupiter.api.Test;

import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Group;
import util.graph.Graph.GraphBuilder;

class GraphTest {
	public static class Label extends AbstractElement<Label> {
		private String label;

		public Label(String label) {
			this.label = label;
		}

		@Override
		public int compareTo(Label o) {
			return label.compareTo(o.label);
		}

		@Override
		public String toString() {
			return label;
		}
	}

	@Test
	void simpleGraphTest() {
		Integers z = Integers.z();
		GraphBuilder<Label, IntE> builder = new GraphBuilder<>(false);
		Label aaron = new Label("Aaron");
		Label alice = new Label("Alice");
		Label bob = new Label("Bob");
		Label charley = new Label("Charley");
		Label delta = new Label("Delta");
		Label eve = new Label("Eve");
		Label mallory = new Label("Mallory");
		Label trent = new Label("Trent");
		Label polly = new Label("Polly");
		Label victor = new Label("Victor");
		Label alpha = new Label("Alpha");
		Label beta = new Label("Beta");
		Label gamma = new Label("Gamma");
		builder.addVertex(aaron);
		builder.addVertex(alice);
		builder.addVertex(bob);
		builder.addVertex(charley);
		builder.addVertex(delta);
		builder.addVertex(eve);
		builder.addVertex(mallory);
		builder.addVertex(trent);
		builder.addVertex(polly);
		builder.addVertex(victor);
		builder.addVertex(alpha);
		builder.addVertex(beta);
		builder.addVertex(gamma);
		builder.addEdge(alice, aaron, z.getInteger(20));
		builder.addEdge(alice, bob, z.getInteger(1));
		builder.addEdge(alice, charley, z.getInteger(3));
		builder.addEdge(charley, bob, z.getInteger(4));
		builder.addEdge(charley, delta, z.getInteger(7));
		builder.addEdge(polly, victor, z.getInteger(2));
		builder.addEdge(trent, victor, z.getInteger(0));
		builder.addEdge(eve, alice, z.getInteger(10));
		builder.addEdge(eve, bob, z.getInteger(10));
		builder.addEdge(eve, charley, z.getInteger(10));
		builder.addEdge(eve, delta, z.getInteger(10));
		builder.addEdge(trent, polly, z.getInteger(5));
		builder.addEdge(trent, delta, z.getInteger(9));
		builder.addEdge(alpha, alpha, z.getInteger(9));
		builder.addEdge(alpha, beta, z.getInteger(3));
		builder.addEdge(gamma, beta, z.getInteger(2));
		builder.addEdge(gamma, beta, z.getInteger(3));
		builder.addEdge(beta, gamma, z.getInteger(1));
		Graph<Label, IntE> graph = builder.build();
		List<Graph<Label, IntE>.Cycle> simple = graph.simpleCycles();
		for (Graph<Label, IntE>.Cycle cycle : simple) {
			System.out.println(cycle);
		}
		List<Graph<Label, IntE>> connected = graph.connectedComponents();
		for (Graph<Label, IntE> component : connected) {
			List<Graph<Label, IntE>.Cycle> simpleCycles = component.simpleCycles();
			if (simpleCycles.isEmpty()) {
				continue;
			}
			Group<Graph<Label, IntE>.Cycle> cycleGroup = component.cycleGroup();
			Graph<Label, IntE>.Cycle twoX = cycleGroup.operate(simpleCycles.get(0), simpleCycles.get(0));
			Graph<Label, IntE>.Cycle XY = cycleGroup.operate(simpleCycles.get(0),
					simpleCycles.get(simpleCycles.size() - 1));
			Graph<Label, IntE>.Cycle X = cycleGroup.operate(twoX, cycleGroup.inverse(simpleCycles.get(0)));
			System.out.println(twoX);
			System.out.println(XY);
			System.out.println(X);
			if (!component.isMultiple()) {
				System.out.println(component.adjacencyMatrix(z.zero()));
			}
		}
		System.out.println(graph.findShortestPath(polly, bob, z.getAdditiveGroup()));
		System.out.println(graph.findShortestPath(alpha, delta, z.getAdditiveGroup()));
		System.out.println(graph.findShortestPath(alpha, gamma, z.getAdditiveGroup()));
	}

}
