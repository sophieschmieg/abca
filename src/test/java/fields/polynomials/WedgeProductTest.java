package fields.polynomials;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;

class WedgeProductTest {

	@Test
	void test() throws IOException {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(f5, 2, Monomial.REVLEX);
		PolynomialRing<PFE> homogenousPolynomials = AbstractPolynomialRing.getPolynomialRing(f5, 3, Monomial.GREVLEX);
		DifferentialForms<PFE> differentialForms = new DifferentialForms<>(polynomials);
		Polynomial<PFE> elliptic = polynomials.parse("-1*Y^2 + X^3 + -1*X");
		System.out.println(differentialForms.derivative(elliptic));
		System.out.println(differentialForms.derivative(differentialForms.derivative(elliptic)));
		Polynomial<PFE> cross = polynomials.parse("X*Y");
		System.out.println(differentialForms.derivative(cross));
		System.out.println(differentialForms.derivative(differentialForms.derivative(cross)));
		ArrayList<Polynomial<PFE>> squareGen = new ArrayList<>();
		squareGen.add(polynomials.parse("X^2"));
		squareGen.add(polynomials.parse("X*Y"));
		squareGen.add(polynomials.parse("Y^2"));
		System.out.println(squareGen);
		PolynomialIdeal<PFE> square = polynomials.getIdeal(squareGen);
		System.out.println(square);
		System.out.println(square.degree());
		List<Polynomial<PFE>> substitute = new ArrayList<>();
		substitute.add(homogenousPolynomials.parse("X + Y + Z"));
		substitute.add(homogenousPolynomials.parse("X + -1*Y"));
		substitute.add(homogenousPolynomials.parse("Z + Y"));
		List<Polynomial<PFE>> newGen = new ArrayList<>();
		for (Polynomial<PFE> gen : square.generators()) {
			newGen.add(polynomials.getEmbedding(homogenousPolynomials
					.dehomogenize(homogenousPolynomials.substitute(polynomials.homogenize(gen), substitute), 3)));
		}
		PolynomialIdeal<PFE> newSquare = polynomials.getIdeal(newGen);
		System.out.println(newSquare);
		System.out.println(newSquare.degree());
		System.out.println(polynomials.primaryDecomposition(square).getPrimaryIdeals());
		System.out.println(polynomials.primaryDecomposition(square).getRadicals());
		System.out.println(polynomials.radical(square));
		System.out.println(polynomials.primaryDecomposition(newSquare).getPrimaryIdeals());
		System.out.println(polynomials.primaryDecomposition(newSquare).getRadicals());
		System.out.println(polynomials.radical(newSquare));
	}

}
