package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.PolynomialRing.GroebnerBasis;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import varieties.affine.AffinePoint;

class PolynomialRingTest {
	
	@Test
	void linearSystemTest() {
		PrimeField fp = PrimeField.getPrimeField(65537);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(fp, 3, Monomial.GREVLEX);
			List<List<Polynomial<PFE>>> equations = new ArrayList<>();
		List<Polynomial<PFE>> row = new ArrayList<>();
		row.add(polynomials.getVar(1));
		row.add(polynomials.getVar(2));
		row.add(polynomials.add(polynomials.getVarPower(1, 3), polynomials.getVar(3)));
		equations.add(row);
		Matrix<Polynomial<PFE>> m = new Matrix<>(equations);
		MatrixModule<Polynomial<PFE>> mm = m.getModule(polynomials);
		System.out.println(mm.kernelBasis(m));
		
	}

	@Test
	void gcdTest() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomials = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
		Polynomial<IntE> p1 = polynomials.subtract(polynomials.getVarPower(1, 2), polynomials.getVarPower(2, 2));
		Polynomial<IntE> p2 = polynomials.add(polynomials.getVarPower(1, 2),
				polynomials.multiply(-2, polynomials.getVar(1), polynomials.getVar(2)), polynomials.getVarPower(2, 2));
		Polynomial<IntE> expected = polynomials.subtract(polynomials.getVar(1), polynomials.getVar(2));
		assertEquals(polynomials.upToUnit(expected), polynomials.upToUnit(polynomials.gcd(p1, p2)));
		Polynomial<IntE> p3 = polynomials.multiply(9, p1);
		Polynomial<IntE> p4 = polynomials.multiply(6, p2);
		Polynomial<IntE> expected2 = polynomials.multiply(3, expected);
		assertEquals(polynomials.upToUnit(expected2), polynomials.upToUnit(polynomials.gcd(p3, p4)));
	}

	@Test
	void solveTest() {
		PrimeField fp = PrimeField.getPrimeField(65537);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(fp, 3, Monomial.GREVLEX);
		Polynomial<PFE> p1 = polynomials.subtract(polynomials.getVarPower(1, 2), polynomials.getVarPower(2, 2));
		Polynomial<PFE> p2 = polynomials.add(polynomials.getVarPower(1, 2),
				polynomials.multiply(-2, polynomials.getVar(1), polynomials.getVar(2)), polynomials.getVarPower(2, 2));
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(p1);
		generators.add(p2);
		generators.add(polynomials.subtract(polynomials.getVarPower(3, 2),
				polynomials.multiply(-12, polynomials.getVarPower(2, 2))));
		generators.add(polynomials.add(polynomials.getVarPower(3, 3), polynomials.getInteger(15)));
		System.out.println(polynomials.solve(generators));
		FiniteField fq = FiniteField.getFiniteField(fp, 6);
		PolynomialRing<FFE> fqPolynomials = AbstractPolynomialRing.getPolynomialRing(fq, 3, Monomial.GREVLEX);
		List<Polynomial<FFE>> embedded = new ArrayList<>();
		for (Polynomial<PFE> p : generators) {
			embedded.add(fqPolynomials.getEmbedding(p, fq.getEmbeddingMap()));
		}
		List<AffinePoint<FFE>> solution = fqPolynomials.solve(embedded);
		for (AffinePoint<FFE> point : solution) {
			System.out.println(point);
			for (Polynomial<FFE> g : embedded) {
				assertEquals(fq.zero(), fqPolynomials.evaluate(g, point));
			}
		}
	}

	@Test
	void buchbergerTest() {
		PrimeField fp = PrimeField.getPrimeField(65537);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(fp, 3, Monomial.GREVLEX);
		Polynomial<PFE> p1 = polynomials.subtract(polynomials.getVarPower(1, 2), polynomials.getVarPower(2, 2));
		Polynomial<PFE> p2 = polynomials.add(polynomials.getVarPower(1, 2),
				polynomials.multiply(-2, polynomials.getVar(1), polynomials.getVar(2)), polynomials.getVarPower(2, 2));
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(p1);
		generators.add(p2);
		Ideal<Polynomial<PFE>> ideal = polynomials.getIdeal(generators);
		List<Polynomial<PFE>> secondSet = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			secondSet.add(ideal.getRandomElement());
		}
		GroebnerBasis<PFE> basis = polynomials.buchberger(secondSet, true);
		for (int i = 0; i < basis.getBasis().size(); i++) {
			Polynomial<PFE> test = polynomials.zero();
			for (int j = 0; j < secondSet.size(); j++) {
				test = polynomials.add(test,
						polynomials.multiply(basis.getExpression().get(i).get(j), secondSet.get(j)));
			}
			assertEquals(basis.getBasis().get(i), test);
		}
	}

	@Test
	void factorizationTest() {
		PrimeField field = PrimeField.getPrimeField(37);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.GREVLEX);
		Polynomial<PFE> p1 = polynomialRing.getLinear(field.one(), field.one());
		Polynomial<PFE> p2 = polynomialRing.add(polynomialRing.getLinear(field.one(), field.negative(field.one())),
				polynomialRing.one());
		Polynomial<PFE> p3 = polynomialRing.add(polynomialRing.multiply(-1, polynomialRing.getVarPower(2, 2)),
				polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1), polynomialRing.one());
		Polynomial<PFE> testPolynomial = polynomialRing.multiply(p1, p2, p3);
		polynomialRing.uniqueFactorization(testPolynomial);
//		for (int tc = 0; tc < 100; tc++) {
//			Polynomial<PFE> p1 = polynomialRing.getRandomElement(3);
//			Polynomial<PFE> p2 = polynomialRing.getRandomElement(2);
//			Polynomial<PFE> p3 = polynomialRing.multiply(p1, p2);
//			polynomialRing.uniqueFactorization(p3);
//			}
	}
}
