package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.PrimaryDecompositionResult;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;

class PrimaryDecompositionTest {

	@Test
	void integersCh13Test() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(1)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(2)));
		generators.add(polynomialRing.getInteger(13));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		assertEquals(3, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(3, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 3; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}

	@Test
	void integersCh0Test() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(1)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(2)));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		assertEquals(3, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(3, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 3; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}

	@Test
	void integersExtraVarTest() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 3, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(1)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(2)));
		generators.add(polynomialRing.getInteger(13));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		System.out.println(ideal);
		for (PolynomialIdeal<IntE> primary : primaryDecomposition.getPrimaryIdeals()) {
			System.out.println(primary);
		}
		assertEquals(3, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(3, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 3; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}
	
	@Test
	void integersDblITestInts() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 2), polynomialRing.one()));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(2, 2), polynomialRing.one()));
		generators.add(polynomialRing.getInteger(30));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		assertEquals(3, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(3, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 3; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}
	
	@Test
	void integersDblITest() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 2), polynomialRing.one()));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(2, 2), polynomialRing.one()));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		assertEquals(2, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(2, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 2; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}
	
	@Test
	void integersDblZIntsTest() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 3, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(3)));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(3)));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		assertEquals(2, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(2, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 2; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}

	@Test
	void integersDblZTest() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 3, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(3)));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(3)));
		generators.add(polynomialRing.getInteger(30));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		PrimaryDecompositionResult<Polynomial<IntE>, PolynomialIdeal<IntE>> primaryDecomposition = polynomialRing
				.primaryDecomposition(ideal);
		assertEquals(3, primaryDecomposition.getPrimaryIdeals().size());
		assertEquals(3, primaryDecomposition.getRadicals().size());
		for (int i = 0; i < 3; i++) {
			assertEquals(primaryDecomposition.getPrimaryIdeals().get(i), primaryDecomposition.getRadicals().get(i));
		}
	}

}
