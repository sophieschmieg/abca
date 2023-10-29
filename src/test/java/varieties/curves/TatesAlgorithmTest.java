package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Field.Extension;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.EllipticCurve.ExtensionLegendreForm;
import varieties.curves.elliptic.EllipticCurve.ExtensionLegendreFormIsomorphism;
import varieties.curves.elliptic.EllipticCurve.Isomorphism;
import varieties.curves.elliptic.EllipticCurve.KodairaSymbol;
import varieties.curves.elliptic.EllipticCurve.TatesAlgorithmResult;
import varieties.projective.ProjectivePoint;

class TatesAlgorithmTest {

	@Test
	void testCoverage() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		int[] a1 = new int[] { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -2 };
		int[] a2 = new int[] { -1, -1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -4, -4, 0, 0,
				32, 8 };
		int[] a3 = new int[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 16, 2 };
		int[] a4 = new int[] { -10, -7820, 0, 4, -36, -171, -1, -2731, -11, -10, -110, 0, 4, -1, -36, -41, -4, -24, -64,
				1, -384, 16, -160, 256, 4, -11, 528, 33 };
		int[] a6 = new int[] { -20, -263580, 0, -6, -70, -874, 0, -55146, 12, -10, -880, 0, 4, 0, -140, -116, 4, -36,
				220, 0, -2772, -180, -1280, -11520, 0, 14, 2432, 38 };
		int[][] primes = new int[][] { { 11 }, { 11 }, { 11 }, { 2, 7 }, { 2, 7 }, { 2, 7 }, { 2, 7 }, { 2, 7 },
				{ 2, 7 }, { 3, 5 }, { 3, 5 }, { 3, 5 }, { 2, 5 }, { 2, 5 }, { 2, 5 }, { 2, 5 }, { 2, 3 }, { 2, 3 },
				{ 2, 3 }, { 2, 3 }, { 2, 3 }, { 2, 3 }, { 2, 11 }, { 2, 3 }, { 2 }, { 2 }, { 2 }, { 2 } };
		KodairaSymbol[][] symbols = new KodairaSymbol[][] { { KodairaSymbol.Iv }, { KodairaSymbol.I1 },
				{ KodairaSymbol.I1 }, { KodairaSymbol.Iv, KodairaSymbol.Iv }, { KodairaSymbol.Iv, KodairaSymbol.Iv },
				{ KodairaSymbol.Iv, KodairaSymbol.I1 }, { KodairaSymbol.I2, KodairaSymbol.I1 },
				{ KodairaSymbol.Iv, KodairaSymbol.I2 }, { KodairaSymbol.I1, KodairaSymbol.I2 },
				{ KodairaSymbol.Iv, KodairaSymbol.Iv }, { KodairaSymbol.Iv, KodairaSymbol.I1 },
				{ KodairaSymbol.I1, KodairaSymbol.I1 }, { KodairaSymbol.IVstar, KodairaSymbol.I2 },
				{ KodairaSymbol.IV, KodairaSymbol.I1 }, { KodairaSymbol.IVstar, KodairaSymbol.Iv },
				{ KodairaSymbol.IV, KodairaSymbol.Iv }, { KodairaSymbol.Ivstar, KodairaSymbol.I2 },
				{ KodairaSymbol.IIIstar, KodairaSymbol.Iv }, { KodairaSymbol.IIIstar, KodairaSymbol.I1 },
				{ KodairaSymbol.III, KodairaSymbol.I1 }, { KodairaSymbol.IIstar, KodairaSymbol.I2 },
				{ KodairaSymbol.IIstar, KodairaSymbol.Iv }, { KodairaSymbol.I0, KodairaSymbol.Iv },
				{ KodairaSymbol.IIstar, KodairaSymbol.Iv }, { KodairaSymbol.Ivstar }, { KodairaSymbol.I0star },
				{ KodairaSymbol.Ivstar }, { KodairaSymbol.Ivstar } };
		int[] conductors = new int[] { 11, 11, 11, 14, 14, 14, 14, 14, 14, 15, 15, 15, 20, 20, 20, 20, 24, 24, 24, 24,
				24, 24, 11, 24, 32, 32, 32, 32 };
		int[][] localValues = new int[][] { { 5 }, { 1 }, { 1 }, { 6, 3 }, { 3, 6 }, { 18, 1 }, { 2, 1 }, { 9, 2 },
				{ 1, 2 }, { 4, 4 }, { 16, 1 }, { 1, 1 }, { -1, 2 }, { -1, 1 }, { -1, 6 }, { -1, 3 }, { 1, 2 },
				{ -1, 4 }, { -1, 1 }, { -1, 1 }, { -1, 2 }, { -1, 8 }, { -1, 5 }, { -1, 8 }, { 3 }, { 0 }, { 3 },
				{ 3 } };
		int[][] localIndeces = new int[][] { { 5 }, { 1 }, { 1 }, { 2, 3 }, { 1, 6 }, { 2, 1 }, { 2, 1 }, { 1, 2 },
				{ 1, 2 }, { 2, 4 }, { 2, 1 }, { 1, 1 }, { 3, 2 }, { 3, 1 }, { 1, 2 }, { 1, 1 }, { 4, 2 }, { 2, 2 },
				{ 2, 1 }, { 2, 1 }, { 1, 2 }, { 1, 2 }, { 1, 5 }, { 1, 2 }, { 4 }, { 2 }, { 4 }, { 4 } };
		boolean[] modelMinimal = new boolean[] { true, true, true, true, true, true, true, true, true, true, true, true,
				true, true, true, true, true, true, true, true, true, true, false, false, true, true, false, true };
		for (int tc = 0; tc < a1.length; tc++) {
			EllipticCurve<Fraction> curve = new EllipticCurve<>(q, q.getInteger(a1[tc]), q.getInteger(a2[tc]),
					q.getInteger(a3[tc]), q.getInteger(a4[tc]), q.getInteger(a6[tc]));
			System.out.println("Testing Curve: " + curve);
			Isomorphism<Fraction> minimalModel = curve.minimalModel(q).get();
			System.out.println("Minimal Model: " + minimalModel.getRange());
			assertEquals(modelMinimal[tc], minimalModel.getScale().equals(q.one()));
			IntE discriminant = curve.discriminant().asInteger();
			IntE conductor = z.one();
			int index = 0;
			for (IntE prime : z.uniqueFactorization(discriminant).primeFactors()) {
				System.out.println("Testing Prime: " + prime);
				assertEquals(z.getInteger(primes[tc][index]), prime);
				TatesAlgorithmResult<Fraction, PFE> reduction = curve.tatesAlgorithm(q.withValuation(prime));
				assertEquals(symbols[tc][index], reduction.getKodairaSymbol());
				if (reduction.getKodairaSymbol().equals(KodairaSymbol.Iv)
						|| reduction.getKodairaSymbol().equals(KodairaSymbol.Ivstar)) {
					assertEquals(localValues[tc][index], reduction.getKodairaSymbolValue());
				}
				assertEquals(localIndeces[tc][index], reduction.getLocalIndex());
				conductor = z.multiply(z.power(prime, reduction.getConductorExponent()), conductor);
				index++;
			}
			assertEquals(conductors[tc], conductor.intValueExact());
		}
	}

	@Test
	void testLegendreFormPrimeField() {
		PrimeField fp = PrimeField.getPrimeField(5);
		for (int tc = 0; tc < 10; tc++) {
			PFE j = fp.getRandomElement();
			PFE scale;
			do {
				scale = fp.getRandomElement();
			} while (scale.equals(fp.zero()));
			EllipticCurve<PFE> curve = new Isomorphism<>(EllipticCurve.fromJInvariant(fp, j), scale,
					fp.getRandomElement(), fp.getRandomElement(), fp.getRandomElement()).getRange();
			System.out.println(curve);
			System.out.println(j);
			ExtensionLegendreForm<PFE, PFE, FFE, FiniteField> legendre = curve
					.getExtensionLengendreForm(fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
			FiniteField ext = legendre.getExtension().extension();
			assertEquals(ext.getEmbedding(curve.jInvariant()), legendre.getCurve().jInvariant());
			for (FFE lambda : legendre.getOtherLambdas()) {
				EllipticCurve<FFE> alternative = EllipticCurve.fromLegendre(ext, lambda);
				assertEquals(ext.getEmbedding(curve.jInvariant()), alternative.jInvariant());

			}
			ExtensionLegendreFormIsomorphism<PFE, PFE, FFE, FiniteField> iso = curve
					.getExtensionLengendreFormIsomorphism(fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
			List<ProjectivePoint<FFE>> torsion = iso.getCurve().getTorsionPoints(2);
			FiniteField field = iso.getExtension().extension();
			for (ProjectivePoint<FFE> point : torsion) {
				if (point.equals(iso.getCurve().neutral())) {
					continue;
				}
				FFE x = point.getDehomogenisedCoord(1, 3);
				assertTrue(x.equals(field.zero()) || x.equals(field.one()) || x.equals(iso.getLambda()));
			}
		}
	}

	//@Test
	void testLegendreFormRationals() {
		Rationals q = Rationals.q();
		for (int tc = 0; tc < 10; tc++) {
			Fraction j = q.getRandomElement();
			j = q.one();
			Fraction scale;
			do {
				scale = q.getRandomElement();
			} while (scale.equals(q.zero()));
			EllipticCurve<Fraction> curve = new Isomorphism<>(EllipticCurve.fromJInvariant(q, j), scale,
					q.getRandomElement(), q.getRandomElement(), q.getRandomElement()).getRange();
			curve = curve.minimalIntegerModel(q).get().getRange();
			System.out.println(curve);
			System.out.println(j);
			ExtensionLegendreForm<Fraction, Fraction, NFE, NumberField> legendre = curve
					.getExtensionLengendreForm(q.getExtension(q.getUnivariatePolynomialRing().getVar()));
			NumberField ext = legendre.getExtension().extension();
			assertEquals(ext.getEmbedding(curve.jInvariant()), legendre.getCurve().jInvariant());
			for (NFE lambda : legendre.getOtherLambdas()) {
				EllipticCurve<NFE> alternative = EllipticCurve.fromLegendre(ext, lambda);
				assertEquals(ext.getEmbedding(curve.jInvariant()), alternative.jInvariant());

			}
			ExtensionLegendreFormIsomorphism<Fraction, Fraction, NFE, NumberField> iso = curve
					.getExtensionLengendreFormIsomorphism(q.getExtension(q.getUnivariatePolynomialRing().getVar()));
			List<ProjectivePoint<NFE>> torsion = iso.getCurve().getTorsionPoints(2);
			NumberField field = iso.getExtension().extension();
			for (ProjectivePoint<NFE> point : torsion) {
				if (point.equals(iso.getCurve().neutral())) {
					continue;
				}
				NFE x = point.getDehomogenisedCoord(1, 3);
				assertTrue(x.equals(field.zero()) || x.equals(field.one()) || x.equals(iso.getLambda()));
			}
		}
	}

	// @Test
	void testStableModel() throws IOException {
		Rationals q = Rationals.q();
		NumberField nf = NumberField.getNumberField(q.getUnivariatePolynomialRing().parse("X^4 + 3"));
		nf = nf.getEmbeddedExtension(nf.getUnivariatePolynomialRing().parse("(X^3 + (-2))")).getField();
		/*
		 * NumberField.getNumberField(q.getUnivariatePolynomialRing().parse("X^2 + -2"))
		 * ; FieldEmbedding<Fraction, NFE, NumberField> embedding = nf
		 * .getEmbeddedExtension(nf.getUnivariatePolynomialRing().parse("X^4 + -3")); nf
		 * = embedding.getField(); NFE sqrt3 = nf.power(embedding.getEmbeddedAlpha(),
		 * 2); UnivariatePolynomialRing<NFE> p = nf.getUnivariatePolynomialRing(); nf =
		 * nf.getEmbeddedExtension( p.toUnivariate(p.subtract(p.getVarPower(2),
		 * p.getEmbedding(nf.add(sqrt3, nf.getInteger(2)))))) .getField();
		 */
		Extension<Fraction, Fraction, NFE, NumberField> ext = q.getExtension(nf.minimalPolynomial());
		EllipticCurve<Fraction> curve1 = new EllipticCurve<>(q, q.zero(), q.getInteger(-1));
		EllipticCurve<NFE> curve1ex = curve1.extendBaseField(ext).getCurve();
		Isomorphism<NFE> isomorphism = new Isomorphism<>(curve1ex, nf.getInteger(2), nf.getInteger(2), nf.zero(),
				nf.getInteger(4));
		System.out.println(isomorphism.getRange());
		System.out.println(curve1.jInvariant());
		System.out.println(curve1.discriminant());
		System.out.println(curve1.getB2());
		System.out.println(curve1.getB4());
		System.out.println(curve1.getB6());
		System.out.println(curve1.getB8());
		System.out.println(curve1.getC4());
		System.out.println(curve1.getC6());
		System.out.println(
				curve1ex.tatesAlgorithm(nf.maximalOrder().localizeAndQuotient(nf.maximalOrder().idealsOver(2).get(0)))
						.getKodairaSymbolString());
		// Y^2 + 8*Y + 16 = X^3 + a*X
		// Y^2 + 8*Y + 16 = (X+r)^3 + a*(X+r)
		// Y^2 + 8*Y + 16 = X^3 + 3rX^2 +3r^2X+r^3 + a*X+r*a
		// Y^2 + 8*Y + 16 = X^3 + 3rX^2 +(3r^2+a)X+ (r^2+a)*r
		//
		// Y^2 + 8*Y + 16 = X^3 + 3rX^2 +(3r^2+a)X+ (r^2+a)*r

		// 3r^2 + a = 16
		// r^3 + a*r = 16

	}

	@Test
	void testHeight() {
		NumberField nf = NumberField.getNumberField();
		EllipticCurve<NFE> curve = new EllipticCurve<>(nf, nf.zero(), nf.getInteger(-2));
		ProjectivePoint<NFE> point = new ProjectivePoint<>(nf, nf.getInteger(3), nf.getInteger(5), nf.one());
		assertTrue(curve.hasRationalPoint(point));
		System.out.println(curve.height(point));
		System.out.println(curve.add(point, point));
		System.out.println(curve.height(curve.add(point, point)));
		System.out.println(curve.neronTatePairing(point, point));
		System.out.println(curve.neronTatePairing(curve.multiply(2,point), point));
			}

}
