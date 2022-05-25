package fields.vectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.ModuloIntegerRing;
import fields.finitefields.ModuloIntegerRing.ModuloIntegerRingElement;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;

public class SmithNormalFormTest<T extends Element<T>> {
	private static <T extends Element<T>> void smithNormalFormTest(Ring<T> ring, int rows, int cols) {
		MatrixModule<T> module = new MatrixModule<>(ring, rows, cols);
		Matrix<T> random = module.getRandomElement();
		System.out.println(random);
		MatrixModule<T>.SmithNormalFormResult decomp = module.smithNormalForm(random);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (i == j) {
					if (i > 0) {
						assertTrue(ring.isDivisible(decomp.getDiagonalMatrix().entry(i + 1, j + 1),
								decomp.getDiagonalMatrix().entry(i, j)));
					}
				} else {
					assertEquals(ring.zero(), decomp.getDiagonalMatrix().entry(i + 1, j + 1));
				}
			}
		}
		System.out.println("R^-1.A.C^-1");
		Matrix<T> rac = module.multiply(module.multiply(decomp.getRowOperationsInverse(), random),
				decomp.getColOperationsInverse());
		System.out.println(rac);
		assertEquals(decomp.getDiagonalMatrix(), rac);
		System.out.println("R.D.C");
		Matrix<T> rdc = module.multiply(module.multiply(decomp.getRowOperations(), decomp.getDiagonalMatrix()),
				decomp.getColOperations());
		System.out.println(rdc);
		assertEquals(random, rdc);
		System.out.println("R.R^-1");
		Matrix<T> rr = module.codomainAlgebra().multiply(decomp.getRowOperations(), decomp.getRowOperationsInverse());
		System.out.println(rr);
		assertEquals(module.codomainAlgebra().one(), rr);
		System.out.println("C^-1.C");
		Matrix<T> cc = module.domainAlgebra().multiply(decomp.getColOperationsInverse(), decomp.getColOperations());
		System.out.println(cc);
		assertEquals(module.domainAlgebra().one(), cc);
		if (rows == cols) {
			MatrixAlgebra<T> algebra = module.domainAlgebra();
			System.out.println("det(A) = " + algebra.determinant(random));
			if (algebra.isUnit(random)) {
				System.out.println("A^-1");
				Matrix<T> inv = algebra.inverse(random);
				System.out.println(inv);
				System.out.println("A^-1.A");
				System.out.println(algebra.multiply(inv, random));
				System.out.println("A.A^-1");
				System.out.println(algebra.multiply(random, inv));
			}
			if (ring.getUnivariatePolynomialRing().isEuclidean()) {
				UnivariatePolynomial<T> ch = algebra.characteristicPolynomial(random);
				System.out.println(ch);
				if (ring instanceof Field) {
					System.out.println(((Field<T>) ring).factorization(ch));
				}
			}
		}
		FreeModule<T> free = module.domain();
		Vector<T> x = free.getRandomElement();
		System.out.println("Ker(A)");
		System.out.println(module.kernelBasis(random));
		System.out.println("x");
		System.out.println(x);
		System.out.println("A.x");
		Vector<T> b = module.multiply(random, x);
		System.out.println(b);
		System.out.println("A.x = b");
		Vector<T> solution = module.solve(random, b);
		System.out.println(solution);
		assertEquals(b, module.multiply(random, solution));
	}

	@Test
	void testPrimeField() {
		PrimeField base = PrimeField.getPrimeField(5);
		UnivariatePolynomialRing<PFE> ring = base.getUnivariatePolynomialRing();
		UnivariatePolynomial<PFE> mini = ring.getPolynomial(base.getElement(2), base.getElement(0), base.getElement(1));
		FiniteField f25 = FiniteField.getFiniteField(mini, base);
		smithNormalFormTest(base, 2, 3);
		smithNormalFormTest(f25, 3, 3);
		smithNormalFormTest(f25.getUnivariatePolynomialRing(), 2, 3);
		smithNormalFormTest(base, 2, 2);
		smithNormalFormTest(f25, 3, 2);
		smithNormalFormTest(f25.getUnivariatePolynomialRing(), 2, 2);
	}

	@Test
	void testRationals() {
		Rationals q = Rationals.q();
		smithNormalFormTest(q, 3, 2);
		smithNormalFormTest(q, 2, 3);
		smithNormalFormTest(q, 3, 3);
		smithNormalFormTest(q, 2, 2);
		smithNormalFormTest(q.getUnivariatePolynomialRing(), 3, 3);
		smithNormalFormTest(q.getUnivariatePolynomialRing(), 3, 2);
		smithNormalFormTest(q.getUnivariatePolynomialRing(), 2, 3);
		smithNormalFormTest(q.getUnivariatePolynomialRing(), 2, 2);
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> mini = ring.getPolynomial(q.getInteger(2), q.getInteger(0), q.getInteger(1));
		NumberField nf = NumberField.getNumberField(mini);
		smithNormalFormTest(nf, 3, 3);
		smithNormalFormTest(nf.getUnivariatePolynomialRing(), 2, 3);
		smithNormalFormTest(nf, 3, 2);
		smithNormalFormTest(nf.getUnivariatePolynomialRing(), 2, 2);
	}

	@Test
	void testIntegers() {
		Integers z = Integers.z();
		List<List<IntE>> m = new ArrayList<>();
		List<IntE> row1 = new ArrayList<>();
		row1.add(z.getInteger(4));
		row1.add(z.zero());
		List<IntE> row2 = new ArrayList<>();
		row2.add(z.zero());
		row2.add(z.getInteger(6));
		m.add(row1);
		m.add(row2);
		Matrix<IntE> matrix = new Matrix<>(m);
		MatrixModule<IntE> mm = matrix.getModule(z);
		mm.smithNormalForm(matrix);
		smithNormalFormTest(Integers.z(), 2, 3);
		smithNormalFormTest(Integers.z(), 3, 3);
		smithNormalFormTest(Integers.z(), 2, 3);
		smithNormalFormTest(Integers.z(), 3, 3);
		smithNormalFormTest(Integers.z(), 2, 2);
		smithNormalFormTest(Integers.z(), 3, 2);
	}

	@Test
	void testIntegersModulo() {
		ModuloIntegerRing mod = new ModuloIntegerRing(625);
		List<List<ModuloIntegerRingElement>> m = new ArrayList<>();
		List<ModuloIntegerRingElement> row1 = new ArrayList<>();
		row1.add(mod.getInteger(25));
		row1.add(mod.getInteger(10));
		List<ModuloIntegerRingElement> row2 = new ArrayList<>();
		row2.add(mod.getInteger(125));
		row2.add(mod.getInteger(25));
		m.add(row1);
		m.add(row2);
		Matrix<ModuloIntegerRingElement> matrix = new Matrix<>(m);
		MatrixModule<ModuloIntegerRingElement> mm = matrix.getModule(mod);
		MatrixModule<ModuloIntegerRingElement>.SmithNormalFormResult smith = mm.smithNormalForm(matrix);
		Matrix<ModuloIntegerRingElement> rac = mm.multiply(mm.multiply(smith.getRowOperationsInverse(), matrix),
				smith.getColOperationsInverse());
		assertEquals(smith.getDiagonalMatrix(), rac);
		Matrix<ModuloIntegerRingElement> rdc = mm
				.multiply(mm.multiply(smith.getRowOperations(), smith.getDiagonalMatrix()), smith.getColOperations());
		assertEquals(matrix, rdc);
		smithNormalFormTest(mod, 2, 3);
		smithNormalFormTest(mod, 3, 3);
		smithNormalFormTest(mod, 2, 3);
		smithNormalFormTest(mod, 3, 3);
		smithNormalFormTest(mod, 2, 2);
		smithNormalFormTest(mod, 3, 2);
	}

	@Test
	void testRealsNotFullRank() {
		// Reals r = Reals.r(128);
		// smithNormalFormTest(r, 2, 3);
	}
}
