package graphics;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.vectors.Matrix;
import fields.vectors.Vector;

class LatticeGraphicsTest {

	@Test
	void test() throws IOException {
		NumberField nf = NumberField.getNumberField(Rationals.q().getUnivariatePolynomialRing().parse("X^2 + 1"));
		Reals r = nf.minkowskiEmbeddingSpace().getValueField();
		Canvas2D<Real, Vector<Real>> canvas = new Canvas2D<>(nf.minkowskiEmbeddingSpace(), 700, 500, r.getDouble(-0.5),
				r.getDouble(6.5), r.getDouble(-0.5), r.getDouble(4.5));
		canvas.setCoordinates(true);
		canvas.setLattice(nf.maximalOrder(), -4, 7, -4, 7);
		NFE point = nf.parse("5/3 + 3/4*i");
		Rationals q = Rationals.q();
		List<Vector<Fraction>> basis = new ArrayList<>();
		basis.add(new Vector<>(q.getInteger(2), q.getInteger(3)));
		basis.add(new Vector<>(q.getInteger(5), q.getInteger(7)));
		Matrix<Fraction> matrix = Matrix.fromColumns(basis);
		System.out.println(matrix.getModule(q).solve(matrix, nf.asVector(point)));
		System.out.println(matrix.getModule(q).solve(matrix, new Vector<>(q.getInteger(2), q.getInteger(1))));
		System.out.println(nf.minkowskiEmbedding(point));
		System.out.println(nf.minkowskiEmbedding(nf.maximalOrder().roundToInteger(point)));
		canvas.addMarkedPoint(nf.minkowskiEmbedding(point), Color.RED);
		canvas.addMarkedPoint(nf.minkowskiEmbedding(nf.parse("4 + 4*i")), Color.PINK);
		canvas.waitUntilClose();
	}

}
