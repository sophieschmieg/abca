package fields.vectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.ModuloIntegerRing;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Ring;
import fields.interfaces.ValueField;

class LinearEquationsTest {
	private static <T extends Element<T>> void equalsVector(Ring<T> ring, FreeModule<T> free, Vector<T> expected,
			Vector<T> actual) {
		if (ring instanceof ValueField<?>) {
			Vector<T> difference = free.subtract(expected, actual);
			ValueField<T> field = (ValueField<T>) ring;
			Reals r = field.getReals();
			Real eps = r.getPowerOfTwo(-r.precision() / 2);
			for (int i = 0; i < difference.dimension(); i++) {
				assertTrue(field.value(difference.get(i + 1)).compareTo(eps) < 0);
			}
		} else {
			assertEquals(expected, actual);
		}
	}

	private static <T extends Element<T>> void linearEquationTest(Ring<T> ring, int rows, int cols, int testcases) {
		MatrixModule<T> module = new MatrixModule<>(ring, rows, cols);
		for (int tc = 0; tc < testcases; tc++) {
			Matrix<T> random = module.getRandomElement();
			System.out.println(random);
			FreeModule<T> free = module.domain();
			Vector<T> x = free.getRandomElement();
			linearEquationTestMatrix(ring, module, random, x);
		}
	}

	private static <T extends Element<T>> void linearEquationTestMatrix(Ring<T> ring, MatrixModule<T> module,
			Matrix<T> matrix, Vector<T> x) {
		System.out.println("Ker(A)");
		List<Vector<T>> kernel = module.kernelBasis(matrix);
		System.out.println(kernel);
		for (Vector<T> kernelVector : kernel) {
			equalsVector(ring, module.codomain(), module.codomain().zero(), module.multiply(matrix, kernelVector));
		}
		Matrix<T> kernelMatrix = null;
		if (kernel.size() > 0) {
			kernelMatrix = Matrix.fromColumns(kernel);
		}
		System.out.println("x");
		System.out.println(x);
		System.out.println("A.x");
		Vector<T> b = module.multiply(matrix, x);
		System.out.println(b);
		System.out.println("A.x = b");
		Vector<T> solution = module.solve(matrix, b);
		System.out.println(solution);
		equalsVector(ring, module.codomain(), b, module.multiply(matrix, solution));
		if (kernelMatrix == null) {
			equalsVector(ring, module.domain(), x, solution);
		} else {
			assertTrue(kernelMatrix.getModule(ring).hasSolution(kernelMatrix, module.domain().subtract(solution, x)));
		}
	}

	@Test
	void testIntegers() {
		Integers z = Integers.z();
		linearEquationTest(z, 2, 3, 100);
		linearEquationTest(z, 3, 2, 100);
		linearEquationTest(z, 2, 2, 100);
		linearEquationTest(z, 3, 3, 100);
	}

	@Test
	void testFractions() {
		Rationals q = Rationals.q();
		List<List<Fraction>> testMatrix = new ArrayList<>();
		List<Fraction> row1 = new ArrayList<>();
		row1.add(q.zero());
		row1.add(q.getInteger(2));
		row1.add(q.getInteger(3));
		row1.add(q.getInteger(4));
		row1.add(q.getInteger(5));
		row1.add(q.getInteger(6));
		testMatrix.add(row1);
		List<Fraction> row2 = new ArrayList<>();
		row2.add(q.zero());
		row2.add(q.zero());
		row2.add(q.zero());
		row2.add(q.getInteger(6));
		row2.add(q.getInteger(5));
		row2.add(q.getInteger(4));
		List<Fraction> row3 = new ArrayList<>();
		row3.add(q.zero());
		row3.add(q.zero());
		row3.add(q.zero());
		row3.add(q.zero());
		row3.add(q.zero());
		row3.add(q.zero());
		testMatrix.add(row3);
		testMatrix.add(row2);
		List<Fraction> testVector = new ArrayList<>();
		testVector.add(q.one());
		testVector.add(q.one());
		testVector.add(q.one());
		testVector.add(q.one());
		testVector.add(q.one());
		testVector.add(q.one());
		linearEquationTestMatrix(q, new MatrixModule<>(q, 3, 6), new Matrix<>(testMatrix), new Vector<>(testVector));
		linearEquationTest(q, 2, 3, 100);
		linearEquationTest(q, 3, 2, 100);
		linearEquationTest(q, 2, 2, 100);
		linearEquationTest(q, 3, 3, 100);
	}

	@Test
	void testReals() {
		Reals r = Reals.r(128);
		linearEquationTest(r, 2, 3, 100);
		linearEquationTest(r, 3, 2, 100);
		linearEquationTest(r, 2, 2, 100);
		linearEquationTest(r, 3, 3, 100);
	}

	@Test
	void testModulo625() {
		ModuloIntegerRing mod = new ModuloIntegerRing(625);
		linearEquationTest(mod, 2, 3, 100);
		linearEquationTest(mod, 3, 2, 100);
		linearEquationTest(mod, 2, 2, 100);
		linearEquationTest(mod, 3, 3, 100);
	}

	@Test
	void testModulo18() {
		ModuloIntegerRing mod = new ModuloIntegerRing(18);
		linearEquationTest(mod, 2, 3, 100);
		linearEquationTest(mod, 3, 2, 100);
		linearEquationTest(mod, 2, 2, 100);
		linearEquationTest(mod, 3, 3, 100);
	}
}
