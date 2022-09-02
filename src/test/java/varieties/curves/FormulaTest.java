package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.QuotientAndRemainderResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.GenericProjectiveScheme;
import varieties.projective.ProjectivePoint;
import varieties.projective.ProjectiveScheme;

class FormulaTest {

	// @Test
	void test() {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q,
				new String[] { "a1", "a2", "a3", "a4", "a6", "s", "x0", "xy0", "y0" });
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

		TranscendentalFieldExtension<Fraction> tr2 = new TranscendentalFieldExtension<>(q,
				new String[] { "a1", "a2", "a3", "a4", "a6" });
		PolynomialRing<TExt<Fraction>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, Monomial.REVLEX,
				new String[] { "s", "x0", "xy0", "y0" });
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

	// @Test
	void testIsomorphisms() {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q,
				AbstractPolynomialRing.getPolynomialRing(q,
						new Monomial.EliminationOrder(Monomial.GREVLEX, Monomial.GREVLEX, 2),
						new String[] { "a1", "a2", "a3", "a4", "a6" }));
		EllipticCurve<TExt<Fraction>> curve = new EllipticCurve<>(tr, tr.getVar(1), tr.getVar(2), tr.getVar(3),
				tr.getVar(4), tr.getVar(5));
		System.out.println(curve);
		System.out.println(curve.getAutomorphisms());
	}

	@Test
	void testDesingularization() throws IOException {
//		Rationals q = Rationals.q();
//		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q,
//				AbstractPolynomialRing.getPolynomialRing(q,
//						new Monomial.EliminationOrder(Monomial.GREVLEX, Monomial.GREVLEX, 2),
//						new String[] { "a0", "a1", "a2", "a3" }));
//		PolynomialRing<TExt<Fraction>> r = AbstractPolynomialRing.getPolynomialRing(tr, Monomial.REVLEX,
//				new String[] { "X", "Z" });
//		Polynomial<TExt<Fraction>> equ = r
//				.parse("((-1)*Z^2 + X^4 + (a3)*X^3*Z + (a2)*X^2*Z^2 + (a1)*X*Z^3 + (a0)*Z^4)");
//		System.out.println(equ);
//		CoordinateRing<TExt<Fraction>> singular = r.getIdeal(Collections.singletonList(equ)).divideOut();
//		CoordinateIdeal<TExt<Fraction>> ideal = (CoordinateIdeal<TExt<Fraction>>) singular.getIdeal(singular.getVar(1),
//				singular.getVar(2));
//		System.out.println(ideal);
//		System.out.println(singular.power(ideal, 2));
//		System.out.println(singular.power(ideal, 3));
//		System.out.println(singular.power(ideal, 4));
//		System.out.println(singular.power(ideal, 5));
		FiniteField fq = FiniteField.getFiniteField(7);
		PolynomialRing<FFE> r = AbstractPolynomialRing.getPolynomialRing(fq, Monomial.REVLEX,
				new String[] { "X", "Y", "V" });
		List<Polynomial<FFE>> list = new ArrayList<>();
		list.add(r.parse("(X^6 + (-1) + (-1)*Y^2)"));
		list.add(r.parse("(V*X + (-1))"));
		PolynomialIdeal<FFE> ideal = r.getIdeal(list);
		CoordinateRing<FFE> cr = ideal.divideOut();
		System.out.println(cr);
		System.out.println(cr.krullDimension());
		AffineScheme<FFE> curve = new AffineScheme<>(fq, cr);
		for (AffinePoint<FFE> point : curve) {
			System.out.println(point);
		}
		System.out.println(curve.singularPoints());
		ProjectiveScheme<FFE> proj = GenericProjectiveScheme.fromAffineScheme(curve);
		for (ProjectivePoint<FFE> point : proj) {
			System.out.println(point);
		}
		System.out.println(proj.singularPoints());
		
		// UnivariatePolynomialRing<FFE> polynomialRing =
		// fq.getUnivariatePolynomialRing();
//		UnivariatePolynomial<FFE> f = polynomialRing
//				.toUnivariate(polynomialRing.subtract(polynomialRing.getVarPower(6), polynomialRing.one()));
//		PolynomialRing<FFE> twoVars = polynomialRing.homogenize(f).getPolynomialRing();
//		UnivariatePolynomial<FFE> fbar = polynomialRing
//				.getEmbedding(twoVars.dehomogenize(polynomialRing.homogenize(f), 1), new int[] { -1, 0 });
//		List<Polynomial<FFE>> polynomials = new ArrayList<>();
//polynomials.add(r.subtract(r.getEmbedding(fq.one(), new int[] {0, 2, 4, 0}), polynomialRing.substitute(fbar, Collections.singletonList(r.getVar(3)))));
////		UnivariatePolynomialRing<FFE> polynomialRing = fq.getUnivariatePolynomialRing();
//		UnivariatePolynomial<FFE> f = polynomialRing.getPolynomial(fq.zero(), fq.getInteger(2), fq.getInteger(-1),
//				fq.getInteger(-2), fq.one());
//		AffineScheme<FFE> firstScheme = new AffineScheme<FFE>(fq,
//				r.getIdeal(Collections.singletonList(r.subtract(r.getVarPower(2, 2), r.getEmbedding(f)))).divideOut());
//		AffineScheme<FFE> secondScheme = new AffineScheme<FFE>(fq, r.getIdeal(Collections.singletonList(
//				r.subtract(r.getVarPower(1, 2), r.dehomogenize(r.getEmbedding(polynomialRing.homogenize(f)), 1))))
//				.divideOut());
//		System.out.println(firstScheme);
//		System.out.println(secondScheme);
//		AffineMorphism<FFE> intersection = AffineScheme.restrictAwayFrom(firstScheme,
//				Collections.singletonList(firstScheme.getCoordinateRing().getVar(1)));
//		CoordinateRing<FFE> intersectionCoordinates = intersection.getDomain().getCoordinateRing();
//		List<CoordinateRingElement<FFE>> mapping = new ArrayList<>();
//		mapping.add(intersectionCoordinates.multiply(intersectionCoordinates.getVar(2),
//				intersectionCoordinates.getVar(3), intersectionCoordinates.getVar(3)));
//		mapping.add(intersectionCoordinates.getVar(3));
//		AffineMorphism<FFE> morphism = new AffineMorphism<>(intersection.getDomain(), secondScheme, mapping);
//		List<AffineScheme<FFE>> cover = new ArrayList<>();
//		cover.add(firstScheme);
//		cover.add(secondScheme);
//		List<List<AffineMorphism<FFE>>> atlas = new ArrayList<>();
//		List<AffineMorphism<FFE>> firstAtlas = new ArrayList<>();
//		firstAtlas.add(firstScheme.identityMorphism());
//		firstAtlas.add(intersection);
//		atlas.add(firstAtlas);
//		List<AffineMorphism<FFE>> secondAtlas = new ArrayList<>();
//		secondAtlas.add(morphism);
//		secondAtlas.add(secondScheme.identityMorphism());
//		atlas.add(secondAtlas);
//		GluedScheme<FFE> glued = new GluedScheme<>(fq, new AffineCover<>(cover, atlas));
//		System.out.println(glued.singularPoints());
	}

	// @Test
	void testSatohSubstitution() {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q,
				AbstractPolynomialRing.getPolynomialRing(q,
						new Monomial.EliminationOrder(Monomial.GREVLEX, Monomial.GREVLEX, 2),
						new String[] { "s", "x0", "a4", "a6", "a4Target", "a6Target" }));
		PolynomialRing<TExt<Fraction>> r = AbstractPolynomialRing.getPolynomialRing(tr, 2, Monomial.LEX);
		Polynomial<TExt<Fraction>> genericRhs = r.add(r.getVarPower(1, 3), r.getVarPower(1, 2),
				r.multiply(tr.getVar(3), r.getVar(1)), r.getEmbedding(tr.getVar(4)));
		Polynomial<TExt<Fraction>> genericLhs = r.getVarPower(2, 2);
		Polynomial<TExt<Fraction>> generic = r.subtract(genericRhs, genericLhs);
		List<Polynomial<TExt<Fraction>>> subs = new ArrayList<>();
		subs.add(r.add(r.multiply(tr.getVarPower(1, -2), r.getVar(1)), r.getEmbedding(tr.getVar(2))));
		subs.add(r.multiply(tr.getVarPower(1, -3), r.getVar(2)));
		Polynomial<TExt<Fraction>> substituted = r.normalize(r.substitute(generic, subs));
		List<Polynomial<Fraction>> coeffs = new ArrayList<>();
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 1, 1 })).asInteger());
		coeffs.add(tr.subtract(substituted.coefficient(r.getMonomial(new int[] { 2, 0 })), tr.one()).asInteger());
		coeffs.add(substituted.coefficient(r.getMonomial(new int[] { 0, 1 })).asInteger());
		coeffs.add(tr.subtract(substituted.coefficient(r.getMonomial(new int[] { 1, 0 })), tr.getVar(5)).asInteger());
		coeffs.add(tr.subtract(substituted.coefficient(r.getMonomial(new int[] { 0, 0 })), tr.getVar(6)).asInteger());
		for (Polynomial<Fraction> coeff : coeffs) {
			System.out.println(coeff);
		}
//		EllipticCurve<TExt<Fraction>> curve1 = new EllipticCurve<>(tr, tr.zero(), tr.one(), tr.zero(), tr.getVar(1),
//				tr.getVar(2));
//		EllipticCurve<TExt<Fraction>> curve2 = new EllipticCurve<>(tr, tr.zero(), tr.one(), tr.zero(), tr.getVar(3),
//				tr.getVar(4));
//		TExt<Fraction> jEquation = tr.subtract(curve1.jInvariant(), curve2.jInvariant());
//		coeffs.add(jEquation.getNumerator());
		PolynomialRing<Fraction> r2 = tr.polynomialRing();
		PolynomialRing<Fraction> r2Elim1 = r2.eliminateVariable(1);
		TranscendentalFieldExtension<Fraction> trSolution1 = new TranscendentalFieldExtension<>(q, r2Elim1);
		PolynomialRing<Fraction> r2Elim2 = r2.eliminateVariable(2);
		TranscendentalFieldExtension<Fraction> trSolution2 = new TranscendentalFieldExtension<>(q, r2Elim2);
		for (Polynomial<Fraction> generator : r2.getIdeal(coeffs).generators()) {
			if (generator.degree(1) == 2) {
				UnivariatePolynomial<Polynomial<Fraction>> regrouped = r2.asUnivariatePolynomial(generator, 1);
				Polynomial<Fraction> denominator = regrouped.univariateCoefficient(2);
				Polynomial<Fraction> assertZero = regrouped.univariateCoefficient(1);
				assertEquals(r2Elim1.zero(), assertZero);
				Polynomial<Fraction> numerator = r2Elim1.negative(regrouped.univariateCoefficient(0));
				TExt<Fraction> solution = trSolution1.getElement(numerator, denominator);
				System.out.println("s^2 = (" + removeDenominators(solution.getNumerator(), r2Elim1) + ")/("
						+ removeDenominators(solution.getDenominator(), r2Elim1) + ")");
			}
			if (generator.degree(1) == 0 && generator.degree(2) == 1) {
				UnivariatePolynomial<Polynomial<Fraction>> regrouped = r2.asUnivariatePolynomial(generator, 2);
				Polynomial<Fraction> denominator = regrouped.univariateCoefficient(1);
				Polynomial<Fraction> numerator = r2Elim2.negative(regrouped.univariateCoefficient(0));
				TExt<Fraction> solution = trSolution2.getElement(numerator, denominator);
				System.out.println("x0 = (" + removeDenominators(solution.getNumerator(), r2Elim2) + ")/("
						+ removeDenominators(solution.getDenominator(), r2Elim2) + ")");
			}
		}
	}

	// @Test
	void testChar2Supersingular() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		TranscendentalFieldExtension<PFE> tr = new TranscendentalFieldExtension<>(f2,
				new String[] { "a2", "a3", "a4", "a6", "s", "x0", "xy0", "y0" });
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

		TranscendentalFieldExtension<PFE> tr2 = new TranscendentalFieldExtension<>(f2,
				new String[] { "a2", "a3", "a4", "a6" });
		PolynomialRing<TExt<PFE>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, Monomial.REVLEX,
				new String[] { "s", "x0", "xy0", "y0" });
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

	// @Test
	void testChar3Supersingular() {
		PrimeField f3 = PrimeField.getPrimeField(3);
		TranscendentalFieldExtension<PFE> tr = new TranscendentalFieldExtension<>(f3,
				new String[] { "a4", "a6", "s", "x0" });
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

		TranscendentalFieldExtension<PFE> tr2 = new TranscendentalFieldExtension<>(f3, new String[] { "a4", "a6" });
		PolynomialRing<TExt<PFE>> tr2Polynomials = AbstractPolynomialRing.getPolynomialRing(tr2, Monomial.REVLEX,
				new String[] { "s", "x0" });
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

	// @Test
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

	// @Test
	void veluTest() throws IOException {
		Rationals q = Rationals.q();
		PolynomialRing<Fraction> polynomialRing = AbstractPolynomialRing.getPolynomialRing(q, Monomial.REVLEX,
				new String[] { "a1", "a2", "a3", "a4", "a6", "xp", "yp", "X", "Y" });
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q, polynomialRing);
		Polynomial<Fraction> syzygy1 = tr.polynomialRing()
				.parse("-1*yp^2 + -1*a1*xp*yp + -1*a3*yp + xp^3 + a2*xp^2 + a4*xp + a6");
		Polynomial<Fraction> syzygy2 = tr.polynomialRing()
				.parse("-1*Y^2 + -1*a1*X*Y + -1*a3*Y + X^3 + a2*X^2 + a4*X + a6");
		List<Polynomial<Fraction>> syzygies = new ArrayList<>();
		syzygies.add(syzygy1);
		syzygies.add(syzygy2);
		PolynomialIdeal<Fraction> ideal = tr.polynomialRing().getIdeal(syzygies);
		TExt<Fraction> s = tr.divide(tr.subtract(tr.getVar(9), tr.getVar(7)), tr.subtract(tr.getVar(8), tr.getVar(6)));
		TExt<Fraction> xr = tr.add(tr.multiply(s, s), tr.negative(tr.getVar(2)), tr.multiply(tr.getVar(1), s),
				tr.negative(tr.add(tr.getVar(6), tr.getVar(8))));
		TExt<Fraction> yThird = tr.add(tr.multiply(s, tr.subtract(xr, tr.getVar(6))), tr.getVar(7));
		TExt<Fraction> yr = tr.negative(tr.add(yThird, tr.multiply(tr.getVar(1), xr), tr.getVar(3)));
		xr = tr.divide(tr.getEmbedding(ideal.residue(xr.getNumerator())),
				tr.getEmbedding(ideal.residue(xr.getDenominator())));
		yThird = tr.divide(tr.getEmbedding(ideal.residue(yThird.getNumerator())),
				tr.getEmbedding(ideal.residue(yThird.getDenominator())));
		yr = tr.divide(tr.getEmbedding(ideal.residue(yr.getNumerator())),
				tr.getEmbedding(ideal.residue(yr.getDenominator())));
		TExt<Fraction> sN = tr.divide(
				tr.add(tr.getVar(9), tr.getVar(7), tr.multiply(tr.getVar(1), tr.getVar(6)), tr.getVar(3)),
				tr.subtract(tr.getVar(8), tr.getVar(6)));
		TExt<Fraction> xrN = tr.add(tr.multiply(sN, sN), tr.negative(tr.getVar(2)), tr.multiply(tr.getVar(1), sN),
				tr.negative(tr.add(tr.getVar(6), tr.getVar(8))));
		TExt<Fraction> yThirdN = tr.add(tr.multiply(sN, tr.subtract(xrN, tr.getVar(6))),
				tr.negative(tr.add(tr.getVar(7), tr.multiply(tr.getVar(1), tr.getVar(6)), tr.getVar(3))));
		TExt<Fraction> yrN = tr.negative(tr.add(yThirdN, tr.multiply(tr.getVar(1), xrN), tr.getVar(3)));
		xrN = tr.divide(tr.getEmbedding(ideal.residue(xrN.getNumerator())),
				tr.getEmbedding(ideal.residue(xrN.getDenominator())));
		yThirdN = tr.divide(tr.getEmbedding(ideal.residue(yThirdN.getNumerator())),
				tr.getEmbedding(ideal.residue(yThirdN.getDenominator())));
		yrN = tr.divide(tr.getEmbedding(ideal.residue(yrN.getNumerator())),
				tr.getEmbedding(ideal.residue(yrN.getDenominator())));
		TExt<Fraction> veluSummandX = tr.add(tr.subtract(xr, tr.getVar(6)), tr.subtract(xrN, tr.getVar(6)));
		TExt<Fraction> veluSummandY = tr.add(tr.subtract(yr, tr.getVar(7)),
				tr.add(yrN, tr.getEmbedding(tr.polynomialRing().parse("yp + a1*xp + a3"))));
		Polynomial<Fraction> xMinusXp = tr.polynomialRing().parse("X + -1*xp");
		QuotientAndRemainderResult<Polynomial<Fraction>> qr = tr.polynomialRing()
				.quotientAndRemainder(veluSummandX.getNumerator(), xMinusXp);
		System.out.println("veluX/(X-xp): " + qr);
		qr = tr.polynomialRing().quotientAndRemainder(veluSummandY.getNumerator(), xMinusXp);
		System.out.println("veluY/(X-xp): " + qr);
		System.out.println("Q veluY/(X-xp): " + qr.getQuotient());
		System.out.println("R veluY/(X-xp): " + qr.getRemainder());
		qr = tr.polynomialRing().quotientAndRemainder(veluSummandY.getNumerator(),
				tr.polynomialRing().power(xMinusXp, 2));
		System.out.println("veluY/(X-xp)^2: " + qr);
		TExt<Fraction> veluX = tr.add(tr.getVar(8), veluSummandX);
		TExt<Fraction> veluY = tr.add(tr.getVar(9), veluSummandY);
		System.out.println(veluX);
		System.out.println(veluY);

		TExt<Fraction> veluSummandXSingle = tr.subtract(xr, tr.getVar(6));
		TExt<Fraction> veluSummandYSingle = tr.subtract(yr, tr.getVar(7));
		qr = tr.polynomialRing().quotientAndRemainder(veluSummandXSingle.getNumerator(), xMinusXp);
		System.out.println("veluXSingle/(X-xp): " + qr);
		qr = tr.polynomialRing().quotientAndRemainder(veluSummandYSingle.getNumerator(), xMinusXp);
		System.out.println("veluYSingle/(X-xp): " + qr);
		qr = tr.polynomialRing().quotientAndRemainder(veluSummandYSingle.getNumerator(),
				tr.polynomialRing().power(xMinusXp, 2));
		System.out.println("veluYSingle/(X-xp)^2: " + qr);
	}

	// @Test
	void testVeluSubstitution() throws IOException {
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q,
				new String[] { "a1", "a2", "a3", "a4", "a6" });
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
		subs.add(r.getVar(1));
		subs.add(r.subtract(r.getVar(2), r.multiply(r.getEmbedding(tr.getEmbedding(q.parse("1/2"))),
				r.add(r.multiply(tr.getVar(1), r.getVar(1)), r.getEmbedding(tr.getVar(3))))));
		System.out.println(subs.get(0));
		System.out.println(subs.get(1));
		Polynomial<TExt<Fraction>> substituted = r.substitute(generic, subs);
		System.out.println(substituted);
		List<Polynomial<TExt<Fraction>>> subs2 = new ArrayList<>();
		subs2.add(r.subtract(r.getVar(1), r.getEmbedding(
				tr.scalarMultiply(q.parse("1/3"), substituted.coefficient(r.getMonomial(new int[] { 2, 0 }))))));
		subs2.add(r.getVar(2));
		System.out.println(subs2.get(0));
		System.out.println(subs2.get(1));
		Polynomial<TExt<Fraction>> substituted2 = r.substitute(substituted, subs2);
		System.out.println(substituted2);
		System.out.println(substituted2.coefficient(r.getMonomial(new int[] { 1, 0 })));
		System.out.println(substituted2.coefficient(r.getMonomial(new int[] { 0, 0 })));
	}

	private Polynomial<Fraction> thirdModularPolynomial(PolynomialRing<Fraction> r) {
		Polynomial<Fraction> phi3 = r.add(r.getVarPower(1, 4), r.getVarPower(2, 4));
		// [1,0] 1855425871872000000000
		// [1,1] -770845966336000000
		// [2,0] 452984832000000
		// [2,1] 8900222976000
		// [2,2] 2587918086
		// [3,0] 36864000
		// [3,1] -1069956
		// [3,2] 2232
		// [3,3] -1
		// [4,0] 1
		phi3 = r.subtract(phi3, r.multiply(r.getVarPower(1, 3), r.getVarPower(2, 3)));
		phi3 = r.add(phi3, r.multiply(2232, r.add(r.multiply(r.getVarPower(1, 3), r.getVarPower(2, 2)),
				r.multiply(r.getVarPower(1, 2), r.getVarPower(2, 3)))));
		phi3 = r.subtract(phi3, r.multiply(1069956,
				r.add(r.multiply(r.getVarPower(1, 3), r.getVar(2)), r.multiply(r.getVar(1), r.getVarPower(2, 3)))));
		phi3 = r.add(phi3, r.multiply(36864000, r.add(r.getVarPower(1, 3), r.getVarPower(2, 3))));
		phi3 = r.add(phi3, r.multiply(BigInteger.valueOf(2587918086L), r.getVarPower(1, 2), r.getVarPower(2, 2)));
		phi3 = r.add(phi3, r.multiply(new BigInteger("8900222976000"),
				r.add(r.multiply(r.getVarPower(1, 2), r.getVar(2)), r.multiply(r.getVar(1), r.getVarPower(2, 2)))));
		phi3 = r.add(phi3,
				r.multiply(new BigInteger("452984832000000"), r.add(r.getVarPower(1, 2), r.getVarPower(2, 2))));
		phi3 = r.subtract(phi3, r.multiply(new BigInteger("770845966336000000"), r.getVar(1), r.getVar(2)));
		phi3 = r.add(phi3, r.multiply(new BigInteger("1855425871872000000000"), r.add(r.getVar(1), r.getVar(2))));
		return phi3;
	}

	private Polynomial<Fraction> removeDenominators(Polynomial<Fraction> t, PolynomialRing<Fraction> r3) {
		Integers z = Integers.z();
		IntE denominator = z.one();
		for (Monomial m : t.monomials()) {
			denominator = z.lcm(t.coefficient(m).getDenominator(), denominator);
		}
		return r3.multiply(denominator, t);
	}

	@SuppressWarnings("unchecked")
	// @Test
	void testSatoh3() throws IOException {
		// char3: Y^2 = X^3 + X^2 - 1/j
		// char2 lifted: Y^2 + XY = X^3 - 36/(j-1728)X - 1/(j-1728)
		// Y = Y' - 1/2X
		// Y'^2 - XY' + 1/4X^2 + XY' - 1/2X^2 = X^3 - 36/(j-1728)X - 1/(j-1728)
		// Y'^2 = X^3 + 1/4X^2 - 36/(j-1728)X - 1/(j-1728)
		// X = 1/4X'
		// Y' = 1/8Y''
		// 1/64 Y''^2 = 1/64X'^3 + 1/64X'^2 - 9/(j-1728)X' - 1/(j-1728)
		// Y''^2 = X'^3 + X'^2 - 576/(j-1728)X' - 64/(j-1728)
		//
		// Model
		// Y^2 = X^3 + X^2 - 576/(j-1728)X - 64/(j-1728)
		Rationals q = Rationals.q();
		TranscendentalFieldExtension<Fraction> tr = new TranscendentalFieldExtension<>(q, new String[] { "j1", "j2" });
		UnivariatePolynomialRing<TExt<Fraction>> r = tr.getUnivariatePolynomialRing();
		TExt<Fraction> coeff1 = tr.inverse(tr.subtract(tr.getVar(1), tr.getInteger(1728)));
		TExt<Fraction> a1 = tr.multiply(-576, coeff1);
		TExt<Fraction> b1 = tr.multiply(-64, coeff1);
		EllipticCurve<TExt<Fraction>> curve1 = new EllipticCurve<>(tr, tr.zero(), tr.one(), tr.zero(), a1, b1);
		UnivariatePolynomial<TExt<Fraction>> divisionPolynomial = r
				.toUnivariate(r.getEmbedding(curve1.getDivisionPolynomial(3)));
		UnivariatePolynomial<TExt<Fraction>> gx = r.getPolynomial(a1, tr.getInteger(2), tr.getInteger(3));
		UnivariatePolynomial<TExt<Fraction>> gy2 = r
				.toUnivariate(r.multiply(4, r.getPolynomial(b1, a1, tr.one(), tr.one())));
		UnivariatePolynomial<TExt<Fraction>> t = r.toUnivariate(r.multiply(2, gx));
		UnivariatePolynomial<TExt<Fraction>> u = r.add(gy2, r.multiply(t, r.getVar()));
		UnivariatePolynomial<TExt<Fraction>> a2 = r.add(r.getEmbedding(a1), r.multiply(-5, t));
		UnivariatePolynomial<TExt<Fraction>> b2 = r
				.toUnivariate(r.add(r.getEmbedding(b1), r.multiply(-7, u), r.multiply(-4, t)));
		TranscendentalFieldExtension<Fraction> tr2 = new TranscendentalFieldExtension<>(q, new String[] { "a", "b" });
		EllipticCurve<TExt<Fraction>> curve2 = new EllipticCurve<>(tr2, tr2.zero(), tr2.one(), tr2.zero(),
				tr2.getVar(1), tr2.getVar(2));
		TExt<Fraction> jInvariant = curve2.jInvariant();
		Polynomial<Fraction> jInvariantNumerator = jInvariant.getNumerator();
		Polynomial<Fraction> jInvariantDenominator = jInvariant.getDenominator();
		PolynomialRing<TExt<Fraction>> r2 = AbstractPolynomialRing.getPolynomialRing(tr, Monomial.GREVLEX,
				new String[] { "a", "b" });
		MathMap<Fraction, TExt<Fraction>> embedding = new MathMap<>() {
			@Override
			public TExt<Fraction> evaluate(Fraction t) {
				return tr.getEmbedding(t);
			}
		};
		List<Polynomial<TExt<Fraction>>> subs = new ArrayList<>();
		subs.add(a2);
		subs.add(b2);
		UnivariatePolynomial<TExt<Fraction>> jInvariantNumeratorActual = r
				.substitute(r2.getEmbedding(jInvariantNumerator, embedding), subs);
		UnivariatePolynomial<TExt<Fraction>> jInvariantDenominatorActual = r
				.substitute(r2.getEmbedding(jInvariantDenominator, embedding), subs);
		UnivariatePolynomial<TExt<Fraction>> jInvariantEquation = r.toUnivariate(
				r.subtract(r.multiply(tr.getVar(2), jInvariantDenominatorActual), jInvariantNumeratorActual));
		System.out.println(jInvariantEquation);
		System.out.println(divisionPolynomial);
		Integers z = Integers.z();
		for (IntE prime : z.setOfPrimes()) {
			System.out.println();
			System.out.println(prime);
			PrimeField fp = PrimeField.getPrimeField(prime);
			TranscendentalFieldExtension<PFE> trP = new TranscendentalFieldExtension<>(fp, new String[] { "j1", "j2" });
			MathMap<Fraction, PFE> reductionMap = new MathMap<>() {

				@Override
				public PFE evaluate(Fraction t) {
					return fp.getFraction(t);
				}
			};
			MathMap<TExt<Fraction>, TExt<PFE>> transReductionMap = new MathMap<>() {

				@Override
				public TExt<PFE> evaluate(TExt<Fraction> t) {
					return trP.getEmbedding(t, reductionMap);
				}
			};
			UnivariatePolynomialRing<TExt<PFE>> rP = trP.getUnivariatePolynomialRing();
			UnivariatePolynomial<TExt<PFE>> divisionPolynomialP = rP.getEmbedding(divisionPolynomial,
					transReductionMap);
			UnivariatePolynomial<TExt<PFE>> jEquationP = rP.getEmbedding(jInvariantEquation, transReductionMap);
			System.out.println(rP.gcd(divisionPolynomialP, jEquationP));
		}
	}
}
