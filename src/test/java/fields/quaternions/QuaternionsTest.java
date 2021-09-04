package fields.quaternions;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.FiniteRationalVectorSpace;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Ring.FactorizationResult;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.quaternions.AbstractQuaternions.Quaternion;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

class QuaternionsTest {

	@Test
	void testInverseAndAssociative() {
		Rationals q = Rationals.q();
		AbstractQuaternions<Fraction> h = RationalQuaternions.quaternions(q.getInteger(-1), q.getInteger(-1));
		for (int i = 0; i < 100; i++) {
			Quaternion<Fraction> t = h.getRandomElement();
			if (t.equals(h.zero())) {
				continue;
			}
			assertEquals(h.one(), h.multiply(t, h.inverse(t)));
			assertEquals(h.one(), h.multiply(h.inverse(t), t));
			Quaternion<Fraction> t2 = h.getRandomElement();
			Quaternion<Fraction> t3 = h.getRandomElement();
			assertEquals(h.multiply(t, h.multiply(t2, t3)), h.multiply(h.multiply(t, t2), t3));
		}
	}

	@Test
	void testConjugateNormTrace() {
		Rationals q = Rationals.q();
		AbstractQuaternions<Fraction> h = RationalQuaternions.quaternions(q.getInteger(-1), q.getInteger(-1));
		for (int i = 0; i < 100; i++) {
			Quaternion<Fraction> t = h.getRandomElement();
			assertEquals(h.getEmbedding(h.reducedNorm(t)), h.multiply(t, h.conjugate(t)));
			assertEquals(h.getEmbedding(h.reducedTrace(t)), h.add(t, h.conjugate(t)));
		}
	}

	@Test
	void testMatrix() {
		Rationals q = Rationals.q();
		AbstractQuaternions<Fraction> h = RationalQuaternions.quaternions(q.getInteger(1), q.getInteger(1));
		MatrixAlgebra<Fraction> m = new FiniteVectorSpace<>(q, 2).matrixAlgebra();
		Matrix<Fraction> o = m.one();
		List<Vector<Fraction>> mi = new ArrayList<>();
		mi.add(new Vector<>(q.one(), q.zero()));
		mi.add(new Vector<>(q.zero(), q.getInteger(-1)));
		Matrix<Fraction> i = Matrix.fromColumns(mi);
		List<Vector<Fraction>> mj = new ArrayList<>();
		mj.add(new Vector<>(q.zero(), q.one()));
		mj.add(new Vector<>(q.getInteger(1), q.zero()));
		Matrix<Fraction> j = Matrix.fromColumns(mj);
		Matrix<Fraction> k = m.multiply(i, j);
		for (int c = 0; c < 100; c++) {
			Quaternion<Fraction> t1 = h.getRandomElement();
			Quaternion<Fraction> t2 = h.getRandomElement();
			Matrix<Fraction> m1 = toMatrix(t1, m, o, i, j, k);
			Matrix<Fraction> m2 = toMatrix(t2, m, o, i, j, k);
			assertEquals(m.multiply(m1, m2), toMatrix(h.multiply(t1, t2), m, o, i, j, k));
			assertEquals(m.multiply(m2, m1), toMatrix(h.multiply(t2, t1), m, o, i, j, k));
		}
	}

	private Matrix<Fraction> toMatrix(Quaternion<Fraction> t, MatrixAlgebra<Fraction> m, Matrix<Fraction> o,
			Matrix<Fraction> i, Matrix<Fraction> j, Matrix<Fraction> k) {
		Matrix<Fraction> result = m.scalarMultiply(t.realPart(), o);
		result = m.add(result, m.scalarMultiply(t.imaginaryPart().get(1), i));
		result = m.add(result, m.scalarMultiply(t.imaginaryPart().get(2), j));
		result = m.add(result, m.scalarMultiply(t.imaginaryPart().get(3), k));
		return result;
	}

	@Test
	public void testSplitRationals() {
		Rationals q = Rationals.q();
		int[][] testCases = new int[][] { { -29, 5 }, { 25, 5 }, { -1, 1 }, { 12, 1 }, { -13, 16 } };
		for (int[] testCase : testCases) {
			RationalQuaternions h = RationalQuaternions.quaternions(q.getInteger(testCase[0]),
					q.getInteger(testCase[1]));
			int[] squares = new int[] { 1, testCase[0], testCase[1], -1 * testCase[0] * testCase[1] };
			assertEquals(q.one(), h.discriminant());
			MatrixAlgebra<Fraction> m = new FiniteRationalVectorSpace(2).matrixAlgebra();
			for (int i = 0; i < 4; i++) {
				Matrix<Fraction> matrix = h.matrixBasis().get(i);
				assertEquals(m.getInteger(squares[i]), m.multiply(matrix, matrix));
			}
			assertEquals(m.zero(), m.add(m.multiply(h.matrixBasis().get(1), h.matrixBasis().get(2)),
					m.multiply(h.matrixBasis().get(2), h.matrixBasis().get(1))));
		}
	}

	@Test
	public void testSplitPAdics() {
		int accuracy = 20;
		int[] primes = new int[] { 2, 3, 5, 7, 11, 13 };
		for (int prime : primes) {
			PAdicField qp = new PAdicField(BigInteger.valueOf(prime), accuracy);
			PAdicField qp2 = qp.withAccuracy(2 * accuracy);
			for (int a = -32; a < 32; a++) {
				for (int b = -32; b < 32; b++) {
		a = -27;
		b = -27;
					if (a == 0 || b == 0) {
						continue;
					}
					PAdicQuaternions h = PAdicQuaternions.quaternions(qp2, qp2.getInteger(a),
							qp2.getInteger(b));
					if (h.isIntegral()) {
						continue;
					}
					System.out.println(h);
					int[] squares = new int[] { 1, a, b, -1 * a * b };
					MatrixAlgebra<PAdicNumber> m = new FiniteVectorSpace<>(qp, 2).matrixAlgebra();
					for (int i = 0; i < 4; i++) {
						Matrix<PAdicNumber> matrix = h.matrixBasis().get(i);
						Matrix<PAdicNumber> square = m.multiply(matrix, matrix);
						Matrix<PAdicNumber> expected = m.getInteger(squares[i]);
						//assertEquals(expected, square);
						assertEquals(qp.round(expected.entry(1, 1), accuracy >> 1),
								qp.round(square.entry(1, 1), accuracy >> 1));
						assertEquals(qp.round(expected.entry(1, 2), accuracy >> 1),
								qp.round(square.entry(1, 2), accuracy >> 1));
						assertEquals(qp.round(expected.entry(2, 1), accuracy >> 1),
								qp.round(square.entry(2, 1), accuracy >> 1));
						assertEquals(qp.round(expected.entry(2, 2), accuracy >> 1),
								qp.round(square.entry(2, 2), accuracy >> 1));
					}
					Matrix<PAdicNumber> abba = m.add(m.multiply(h.matrixBasis().get(1), h.matrixBasis().get(2)),
							m.multiply(h.matrixBasis().get(2), h.matrixBasis().get(1)));
					//assertEquals(m.zero(), abba);
					assertEquals(qp.zero(), qp.round(abba.entry(1, 1), accuracy >> 1));
					assertEquals(qp.zero(), qp.round(abba.entry(1, 2), accuracy >> 1));
					assertEquals(qp.zero(), qp.round(abba.entry(2, 1), accuracy >> 1));
					assertEquals(qp.zero(), qp.round(abba.entry(2, 2), accuracy >> 1));
				}
			}
		}
	}

	@Test
	public void testDiscriminant() {
		Integers z = Integers.z();
		for (int i = 1; i < 100; i++) {
			IntE iInteger = z.getInteger(i);
			FactorizationResult<IntE> factors = z.uniqueFactorization(iInteger);
			boolean squareFree = true;
			for (IntE prime : factors.primeFactors()) {
				if (factors.multiplicity(prime) > 1) {
					squareFree = false;
					break;
				}
			}
			if (!squareFree) {
				continue;
			}
			RationalQuaternions h = RationalQuaternions.quaternions(iInteger);
			assertEquals(iInteger, h.discriminant().asInteger());

		}
	}
}
