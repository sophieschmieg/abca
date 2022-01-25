package fields.vectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.integers.Integers;
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
		if (!ring.isEuclidean()) {
			return;
		}
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
		NumberField nf = new NumberField(mini);
		smithNormalFormTest(nf, 3, 3);
		smithNormalFormTest(nf.getUnivariatePolynomialRing(), 2, 3);
		smithNormalFormTest(nf, 3, 2);
		smithNormalFormTest(nf.getUnivariatePolynomialRing(), 2, 2);
	}

	@Test
	void testIntegers() {
		smithNormalFormTest(Integers.z(), 2, 3);
		smithNormalFormTest(Integers.z(), 3, 3);
		smithNormalFormTest(Integers.z(), 2, 3);
		smithNormalFormTest(Integers.z(), 3, 3);
		smithNormalFormTest(Integers.z(), 2, 2);
		smithNormalFormTest(Integers.z(), 3, 2);
	}

	@Test
	void testRealsNotFullRank() {
		//Reals r = Reals.r(128);
		//smithNormalFormTest(r, 2, 3);
	}
}
