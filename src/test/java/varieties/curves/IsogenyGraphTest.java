package varieties.curves;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import util.graph.Edge;
import util.graph.Graph;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.Isogeny;

class IsogenyGraphTest {

	@Test
	void isogenyGraphTest() {
		FiniteField fq = FiniteField.getFiniteField(43 * 43);
		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(1728));
		Graph<EllipticCurve<FFE>, Isogeny<FFE>> graph = curve.isogenyGraph(2);
		for (EllipticCurve<FFE> iso : graph.vertexSet()) {
			for (Edge<EllipticCurve<FFE>, Isogeny<FFE>> edge : graph.edges(iso)) {
				System.out.println(edge);
			}
		}
		System.out.println("connected: " + graph.isConnected());
		List<Graph<EllipticCurve<FFE>, Isogeny<FFE>>.Cycle> cycles = graph.simpleCycles();
		for (Graph<EllipticCurve<FFE>, Isogeny<FFE>>.Cycle cycle : cycles) {
			List<String> asStrings= new ArrayList<>();
			for (EllipticCurve<FFE> t : cycle.vertexPath()) {
				asStrings.add(t.jInvariant().toString());
			}
			System.out.println(String.join(", ", asStrings));
		}
	}

}
