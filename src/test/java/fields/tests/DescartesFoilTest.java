package fields.tests;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import varieties.SimpleFunctionField;
import varieties.SimpleFunctionField.SimpleRationalFunction;

class DescartesFoilTest {

	@Test
	void test() throws IOException {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> field = new TranscendentalFieldExtension<>(q, new String[] { "x0" });
		UnivariatePolynomialRing<TExt<Fraction>> univariate = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<TExt<Fraction>> minimalPolynomial = univariate.getPolynomial(field.getVarPower(1, 3),
				field.multiply(-3, field.getVar(1)), field.zero(), field.one());
		SimpleFunctionField<Fraction> functionField = new SimpleFunctionField<>(minimalPolynomial, field, "y0");
		PolynomialRing<SimpleRationalFunction<Fraction>> polynomials = AbstractPolynomialRing
				.getPolynomialRing(functionField, Monomial.GREVLEX, new String[] { "X", "Y", "U" });
		List<Polynomial<SimpleRationalFunction<Fraction>>> generators = new ArrayList<>();
		generators.add(polynomials.parse("X^3 + Y^3 + (-3)*X*Y"));
		Polynomial<SimpleRationalFunction<Fraction>> circle = polynomials.parse("X^2 + (-2)*X*U + Y^2");
		circle = polynomials.add(circle, polynomials.multiply(2,
				polynomials.getEmbedding(functionField.getEmbedding(field.getVar(1))), polynomials.getVar(3)));
		circle = polynomials.subtract(circle,
				polynomials.getEmbedding(functionField.add(functionField.getEmbedding(field.getVarPower(1, 2)),
						functionField.fromPolynomial(field.getUnivariatePolynomialRing().getVarPower(2)))));
		generators.add(circle);
		// generators.add(polynomials.parse("x0^2 + -2*x0*U + U^2 + y0^2 + -1*R"));

		// generators.add(polynomials.parse("Y + -1*y0"));
		System.out.println(generators);
		PolynomialIdeal<SimpleRationalFunction<Fraction>> ideal = polynomials.getIdeal(generators);
		System.out.println(ideal);
		System.out.println(ideal.dimension());
		System.out.println(ideal.degree());
		CoordinateRing<SimpleRationalFunction<Fraction>> cr = ideal.divideOut();
		CoordinateIdeal<SimpleRationalFunction<Fraction>> coord = cr
				.getIdeal(Collections.singletonList(cr.getEmbedding((polynomials.subtract(polynomials.getVar(1),
						polynomials.getEmbedding(functionField.getEmbedding(field.getVar(1))))))));
		LocalizedCoordinateRing<SimpleRationalFunction<Fraction>> local = new LocalizedCoordinateRing<>(functionField, cr, coord);
	System.out.println(local);
	}

}
