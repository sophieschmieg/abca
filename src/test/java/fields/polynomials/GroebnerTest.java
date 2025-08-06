package fields.polynomials;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;

class GroebnerTest {

	@Test
	void test() throws IOException {
		Rationals q = Rationals.q();
		PolynomialRing<Fraction> ring = AbstractPolynomialRing.getPolynomialRing(q, Monomial.GREVLEX,
				new String[] { "X", "Y", "Z" });
		List<Polynomial<Fraction>> list = new ArrayList<>();
		list.add(ring.parse("X*Y + Z + -40"));
		list.add(ring.parse("X*Z + Y + -51"));
		list.add(ring.parse("X + Y + Z + -19"));
		System.out.println(list);
		PolynomialIdeal<Fraction> ideal = ring.getIdeal(list);
		System.out.println(ideal);
		System.out.println(ring.solve(ideal));
	}

}
