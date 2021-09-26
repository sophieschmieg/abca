package fields.helper;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;

class CoordinateRingTest {

	@Test
	void testDegree3() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(1)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		assertEquals(3, coordinateRing.degree());
	}

	@Test
	void testDegreeElliptic() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1), polynomialRing.one()), polynomialRing.getVarPower(2, 2)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		assertEquals(5, coordinateRing.degree());
	}

	@Test
	void testDegreeCommon() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 4, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1)), polynomialRing.getVarPower(2, 2)));
		generators.add(polynomialRing.subtract(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(4)),polynomialRing.multiply(polynomialRing.getVar(2), polynomialRing.getVar(3))));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 4), polynomialRing.getVarPower(2, 4), polynomialRing.getVarPower(3, 2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		int degree = coordinateRing.degree();
		generators.add(polynomialRing.getVar(4));
		PolynomialIdeal<PFE> ideal2 = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing2 = new CoordinateRing<>(polynomialRing, ideal2);
		assertEquals(degree, coordinateRing2.degree());
	}

}
