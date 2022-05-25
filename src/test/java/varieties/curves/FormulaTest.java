package varieties.curves;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;

class FormulaTest {

	@Test
	void test() {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q, 9);
		PolynomialRing<TExt<Fraction>> r = AbstractPolynomialRing.getPolynomialRing(tr, 2, Monomial.LEX);
		Polynomial<TExt<Fraction>> genericRhs = r.add(r.getVarPower(1, 3),
				r.multiply(tr.getVar(2), r.getVarPower(1, 2)), r.multiply(tr.getVar(4), r.getVar(1)),
				r.getEmbedding(tr.getVar(5)));
		Polynomial<TExt<Fraction>> genericLhs = r.add(r.getVarPower(2, 2),
				r.multiply(r.getEmbedding(tr.getVar(1)), r.getVar(1), r.getVar(2)),
				r.multiply(tr.getVar(3), r.getVar(2)));
		Polynomial<TExt<Fraction>> generic = r.subtract(genericRhs, genericLhs);
		System.out.println(generic);
		List<Polynomial<TExt<Fraction>>> subs = new ArrayList<>();
		subs.add(r.add(r.multiply(tr.getVarPower(6, -2), r.getVar(1)), r.getEmbedding(tr.getVar(7))));
		subs.add(r.add(r.multiply(tr.getVarPower(6, -3), r.getVar(2)), r.multiply(tr.getVar(8), r.getVar(1)),
				r.getEmbedding(tr.getVar(9))));
		System.out.println(subs.get(0));
		System.out.println(subs.get(1));
		Polynomial<TExt<Fraction>> substituted = r.normalize(r.substitute(generic, subs));
		System.out.println(substituted);
		System.out.println("a1 = " + tr.negative(substituted.coefficient(r.getMonomial(new int[] { 1, 1 }))));
		System.out.println("a2 = " + substituted.coefficient(r.getMonomial(new int[] { 2, 0 })));
		System.out.println("a3 = " + tr.negative(substituted.coefficient(r.getMonomial(new int[] { 0, 1 }))));
		System.out.println("a4 = " + substituted.coefficient(r.getMonomial(new int[] { 1, 0 })));
		System.out.println("a6 = " + substituted.coefficient(r.getMonomial(new int[] { 0, 0 })));
	}

}
