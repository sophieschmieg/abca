package varieties.curves;

import static org.junit.jupiter.api.Assertions.fail;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;

class FormulaTest {

	@Test
	void test() {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q, new String[] {"a1", "a2", "a3", "a4", "a6", "s", "x0", "xy0", "y0"});
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
		Map<Integer, Polynomial<Fraction>> coeffs = new TreeMap<>();
		coeffs.put(1, tr.negative(substituted.coefficient(r.getMonomial(new int[] { 1, 1 }))).asInteger());
		coeffs.put(2, substituted.coefficient(r.getMonomial(new int[] { 2, 0 })).asInteger());
		coeffs.put(3, tr.negative(substituted.coefficient(r.getMonomial(new int[] { 0, 1 }))).asInteger());
		coeffs.put(4, substituted.coefficient(r.getMonomial(new int[] { 1, 0 })).asInteger());
		coeffs.put(5, substituted.coefficient(r.getMonomial(new int[] { 0, 0 })).asInteger());

		for (int i : coeffs.keySet()) {
			System.out.println("a" + i + " = " + coeffs.get(i));
		}

		TranscendentalFieldExtension<Fraction> tr2 = new TranscendentalFieldExtension<>(q, new String[] {"a1", "a2", "a3", "a4", "a6"});
		PolynomialRing<TExt<Fraction>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, 
				Monomial.REVLEX, new String[] {"s", "x0", "xy0", "y0"});
		List<Polynomial<TExt<Fraction>>> asEquations = new ArrayList<>();
		for (int i : coeffs.keySet()) {
			Polynomial<Fraction> p = coeffs.get(i);
			Polynomial<TExt<Fraction>> transformed = tr2Polynomials
					.negative(tr2Polynomials.getEmbedding(tr2.getVar(i)));
			for (Monomial m : p.monomials()) {
				Polynomial<TExt<Fraction>> coeff = tr2Polynomials.getEmbedding(tr2.getEmbedding(p.coefficient(m)));
				for (int j = 0; j < 5; j++) {
					coeff = tr2Polynomials.multiply(tr2.getVarPower(j + 1, m.exponents()[j]), coeff);
				}
				for (int j = 0; j < 4; j++) {
					coeff = tr2Polynomials.multiply(tr2Polynomials.getVarPower(j + 1, m.exponents()[j + 5]), coeff);
				}
				transformed = tr2Polynomials.add(coeff, transformed);
			}
			asEquations.add(transformed);
		}
		fail();
		System.out.println(tr2Polynomials.getIdeal(asEquations));
	}

	@Test
	void testChar2Supersingular() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		TranscendentalFieldExtension<PFE> tr = new TranscendentalFieldExtension<>(f2, new String[] {"a2", "a3", "a4", "a6", "s", "x0", "xy0", "y0"});
		PolynomialRing<TExt<PFE>> r = AbstractPolynomialRing.getPolynomialRing(tr, 2, Monomial.LEX);
		Polynomial<TExt<PFE>> genericRhs = r.add(r.getVarPower(1, 3), r.multiply(tr.getVar(1), r.getVarPower(1, 2)),
				r.multiply(tr.getVar(3), r.getVar(1)), r.getEmbedding(tr.getVar(4)));
		Polynomial<TExt<PFE>> genericLhs = r.add(r.getVarPower(2, 2), r.multiply(tr.getVar(2), r.getVar(2)));
		Polynomial<TExt<PFE>> generic = r.subtract(genericRhs, genericLhs);
		System.out.println(generic);
		List<Polynomial<TExt<PFE>>> subs = new ArrayList<>();
		subs.add(r.add(r.multiply(tr.getVarPower(5, -2), r.getVar(1)), r.getEmbedding(tr.getVar(6))));
		subs.add(r.add(r.multiply(tr.getVarPower(5, -3), r.getVar(2)), r.multiply(tr.getVar(7), r.getVar(1)),
				r.getEmbedding(tr.getVar(8))));
		System.out.println(subs.get(0));
		System.out.println(subs.get(1));
		Polynomial<TExt<PFE>> substituted = r.normalize(r.substitute(generic, subs));
		System.out.println(substituted);
		List<Polynomial<PFE>> coeffs = new ArrayList<>();
		coeffs.add(tr.negative(substituted.coefficient(r.getMonomial(new int[] { 1, 1 }))).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 2, 0 })).asInteger());
		coeffs.add(tr.negative(substituted.coefficient(r.getMonomial(new int[] { 0, 1 }))).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 1, 0 })).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 0, 0 })).asInteger());

		for (int i = 0; i < coeffs.size(); i++) {
			System.out.println("a" + i + " = " + coeffs.get(i));
		}

		TranscendentalFieldExtension<PFE> tr2 = new TranscendentalFieldExtension<>(f2, new String[] {"a2", "a3", "a4", "a6"});
		PolynomialRing<TExt<PFE>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, Monomial.REVLEX, new String[] {"s", "x0", "xy0", "y0"});
		List<Polynomial<TExt<PFE>>> asEquations = new ArrayList<>();
		for (Polynomial<PFE> p : coeffs) {
			Polynomial<TExt<PFE>> transformed = tr2Polynomials.zero();
			for (Monomial m : p.monomials()) {
				Polynomial<TExt<PFE>> coeff = tr2Polynomials.getEmbedding(tr2.getEmbedding(p.coefficient(m)));
				for (int j = 0; j < 4; j++) {
					coeff = tr2Polynomials.multiply(tr2.getVarPower(j + 1, m.exponents()[j]), coeff);
				}
				for (int j = 0; j < 4; j++) {
					coeff = tr2Polynomials.multiply(tr2Polynomials.getVarPower(j + 1, m.exponents()[j + 4]), coeff);
				}
				transformed = tr2Polynomials.add(coeff, transformed);
			}
			asEquations.add(transformed);
		}
		for (Polynomial<TExt<PFE>> equation : asEquations) {
			System.out.println(equation);
		}
		List<Polynomial<TExt<PFE>>> toSolve = new ArrayList<>();
		toSolve.add(asEquations.get(0));
		toSolve.add(asEquations.get(1));
		toSolve.add(tr2Polynomials.add(asEquations.get(2), tr2Polynomials.one()));
		toSolve.add(asEquations.get(3));
		toSolve.add(asEquations.get(4));
		PolynomialIdeal<TExt<PFE>> ideal = tr2Polynomials.getIdeal(toSolve);
		System.out.println(ideal);
		System.out.println(ideal.dimension());
		System.out.println(ideal.degree());
	}

	@Test
	void testChar3Supersingular() {
		PrimeField f3 = PrimeField.getPrimeField(3);
		TranscendentalFieldExtension<PFE> tr = new TranscendentalFieldExtension<>(f3, new String[] {"a4", "a6", "s", "x0"});
		PolynomialRing<TExt<PFE>> r = AbstractPolynomialRing.getPolynomialRing(tr, 2, Monomial.LEX);
		Polynomial<TExt<PFE>> genericRhs = r.add(r.getVarPower(1, 3), r.multiply(tr.getVar(1), r.getVar(1)),
				r.getEmbedding(tr.getVar(2)));
		Polynomial<TExt<PFE>> genericLhs = r.getVarPower(2, 2);
		Polynomial<TExt<PFE>> generic = r.subtract(genericRhs, genericLhs);
		System.out.println(generic);
		List<Polynomial<TExt<PFE>>> subs = new ArrayList<>();
		subs.add(r.add(r.multiply(tr.getVarPower(3, -2), r.getVar(1)), r.getEmbedding(tr.getVar(4))));
		subs.add(r.multiply(tr.getVarPower(3, -3), r.getVar(2)));
		System.out.println(subs.get(0));
		System.out.println(subs.get(1));
		Polynomial<TExt<PFE>> substituted = r.normalize(r.substitute(generic, subs));
		System.out.println(substituted);
		List<Polynomial<PFE>> coeffs = new ArrayList<>();
		coeffs.add(tr.negative(substituted.coefficient(r.getMonomial(new int[] { 1, 1 }))).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 2, 0 })).asInteger());
		coeffs.add(tr.negative(substituted.coefficient(r.getMonomial(new int[] { 0, 1 }))).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 1, 0 })).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 0, 0 })).asInteger());

		for (int i = 0; i < coeffs.size(); i++) {
			System.out.println("a" + i + " = " + coeffs.get(i));
		}

		TranscendentalFieldExtension<PFE> tr2 = new TranscendentalFieldExtension<>(f3, new String[] {"a4", "a6"});
		PolynomialRing<TExt<PFE>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, Monomial.REVLEX, new String[] {"s", "x0"});
		List<Polynomial<TExt<PFE>>> asEquations = new ArrayList<>();
		for (Polynomial<PFE> p : coeffs) {
			Polynomial<TExt<PFE>> transformed = tr2Polynomials.zero();
			for (Monomial m : p.monomials()) {
				Polynomial<TExt<PFE>> coeff = tr2Polynomials.getEmbedding(tr2.getEmbedding(p.coefficient(m)));
				for (int j = 0; j < 2; j++) {
					coeff = tr2Polynomials.multiply(tr2.getVarPower(j + 1, m.exponents()[j]), coeff);
				}
				for (int j = 0; j < 2; j++) {
					coeff = tr2Polynomials.multiply(tr2Polynomials.getVarPower(j + 1, m.exponents()[j + 2]), coeff);
				}
				transformed = tr2Polynomials.add(coeff, transformed);
			}
			asEquations.add(transformed);
		}
		for (Polynomial<TExt<PFE>> equation : asEquations) {
			System.out.println(equation);
		}
		List<Polynomial<TExt<PFE>>> toSolve = new ArrayList<>();
		toSolve.add(asEquations.get(0));
		toSolve.add(asEquations.get(1));
		toSolve.add(asEquations.get(2));
		toSolve.add(tr2Polynomials.subtract(asEquations.get(3), tr2Polynomials.one()));
		toSolve.add(asEquations.get(4));
		PolynomialIdeal<TExt<PFE>> ideal = tr2Polynomials.getIdeal(toSolve);
		System.out.println(ideal);
		System.out.println(ideal.dimension());
		System.out.println(ideal.degree());
	}

	@Test
	void testChar3() {
		PrimeField f3 = PrimeField.getPrimeField(3);
		TranscendentalFieldExtension<PFE> tr = new TranscendentalFieldExtension<>(f3,
				new String[] { "a2", "a4", "a6", "s", "x0" });
		PolynomialRing<TExt<PFE>> r = AbstractPolynomialRing.getPolynomialRing(tr, 2, Monomial.LEX);
		Polynomial<TExt<PFE>> genericRhs = r.add(r.getVarPower(1, 3), r.multiply(tr.getVar(1), r.getVarPower(1, 2)),
				r.multiply(tr.getVar(2), r.getVar(1)), r.getEmbedding(tr.getVar(3)));
		Polynomial<TExt<PFE>> genericLhs = r.getVarPower(2, 2);
		Polynomial<TExt<PFE>> generic = r.subtract(genericRhs, genericLhs);
		System.out.println(generic);
		List<Polynomial<TExt<PFE>>> subs = new ArrayList<>();
		subs.add(r.add(r.multiply(tr.getVarPower(4, -2), r.getVar(1)), r.getEmbedding(tr.getVar(5))));
		subs.add(r.multiply(tr.getVarPower(4, -3), r.getVar(2)));
		System.out.println(subs.get(0));
		System.out.println(subs.get(1));
		Polynomial<TExt<PFE>> substituted = r.normalize(r.substitute(generic, subs));
		System.out.println(substituted);
		List<Polynomial<PFE>> coeffs = new ArrayList<>();
		coeffs.add(tr.negative(substituted.coefficient(r.getMonomial(new int[] { 1, 1 }))).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 2, 0 })).asInteger());
		coeffs.add(tr.negative(substituted.coefficient(r.getMonomial(new int[] { 0, 1 }))).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 1, 0 })).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 0, 0 })).asInteger());

		int[] indeces = new int[] { 1, 2, 3, 4, 6 };
		for (int i = 0; i < coeffs.size(); i++) {
			System.out.println("a" + indeces[i] + "' = " + coeffs.get(i));
		}

		TranscendentalFieldExtension<PFE> tr2 = new TranscendentalFieldExtension<>(f3,
				new String[] { "a2", "a4", "a6" });
		PolynomialRing<TExt<PFE>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, Monomial.REVLEX,
				new String[] { "s", "x0" });
		List<Polynomial<TExt<PFE>>> asEquations = new ArrayList<>();
		for (Polynomial<PFE> p : coeffs) {
			Polynomial<TExt<PFE>> transformed = tr2Polynomials.zero();
			for (Monomial m : p.monomials()) {
				Polynomial<TExt<PFE>> coeff = tr2Polynomials.getEmbedding(tr2.getEmbedding(p.coefficient(m)));
				for (int j = 0; j < 3; j++) {
					coeff = tr2Polynomials.multiply(tr2.getVarPower(j + 1, m.exponents()[j]), coeff);
				}
				for (int j = 0; j < 2; j++) {
					coeff = tr2Polynomials.multiply(tr2Polynomials.getVarPower(j + 1, m.exponents()[j + 3]), coeff);
				}
				transformed = tr2Polynomials.add(coeff, transformed);
			}
			asEquations.add(transformed);
		}
		for (Polynomial<TExt<PFE>> equation : asEquations) {
			System.out.println(equation);
		}
		List<Polynomial<TExt<PFE>>> toSolve = new ArrayList<>();
		toSolve.add(asEquations.get(0));
		toSolve.add(tr2Polynomials.subtract(asEquations.get(1), tr2Polynomials.one()));
		toSolve.add(asEquations.get(2));
		toSolve.add(asEquations.get(3));
		toSolve.add(asEquations.get(4));
		PolynomialIdeal<TExt<PFE>> ideal = tr2Polynomials.getIdeal(toSolve);
		System.out.println(ideal);
		System.out.println(ideal.dimension());
		System.out.println(ideal.degree());
	}

}
