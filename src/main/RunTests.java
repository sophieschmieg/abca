package main;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PrimeFieldElement;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.helper.ExtensionField;
import fields.helper.ExtensionField.ExtensionFieldElement;
import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.padics.PAdicExtensionField.PAdicNumberExt;
import fields.padics.PAdicField;
import fields.padics.PAdicField.PAdicNumber;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.UnivariatePolynomial;
import fields.polynomials.UnivariatePolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.MiscAlgorithms;
import varieties.AffinePoint;
import varieties.FunctionField;
import varieties.ProjectivePoint;
import varieties.RationalFunction;
import varieties.curves.EllipticCurve;

public class RunTests {

	public static void main(String[] args) {
		PAdicField z5 = new PAdicField(BigInteger.valueOf(5), 12);
		for (int i = 0; i < 100; i++) {
			for (int j = 0; j < 100; j++) {
				PAdicNumber n = z5.getInteger(i);
				PAdicNumber m = z5.getInteger(j);
				if (z5.toRational(n).getNumerator().getValue().intValueExact() != i) {
					System.out.println(n + " != " + i);
				}
				if (z5.toRational(m).getNumerator().getValue().intValueExact() != j) {
					System.out.println(n + " != " + i);
				}
				if (z5.toRational(z5.add(n, m)).getNumerator().getValue().intValueExact() != i + j) {
					System.out.println(n + " + " + m + " != " + i + " + " + j);
				}
				if (z5.toRational(z5.multiply(n, m)).getNumerator().getValue().intValueExact() != i * j) {
					System.out.println(n + " * " + m + " = " + z5.multiply(n, m) + " != " + i + " * " + j);
				}
				if (j != 0 && !z5.multiply(m, z5.divide(n, m)).equals(n)) {
					System.out.println(n + " * (" + m + "/" + n + ") = " + n + " * " + z5.divide(m, n) + " = " + z5.multiply(m, z5.divide(n, m)) + " != " + n);
				}
			}
		}
		System.out.println();
		Rationals q = Rationals.q();
		PolynomialRing<Fraction> q4 = AbstractPolynomialRing.getPolynomialRing(q, 4, Monomial.LEX);
		Polynomial<Fraction> p1 = q4.getVarPower(1, 2);
		p1 = q4.subtract(p1, q4.getVarPower(2, 2));
		p1 = q4.subtract(p1, q4.getVarPower(3, 2));
		p1 = q4.add(p1, q4.getVarPower(4, 2));
		p1 = q4.subtract(p1, q4.getInteger(-1));
		Polynomial<Fraction> p2 = q4.negative(q4.getVarPower(2, 2));
		p2 = q4.add(p2, q4.getVarPower(4, 2));
		p2 = q4.add(p2, q4.multiply(2, q4.getVar(1), q4.getVar(2)));
		p2 = q4.subtract(p2, q4.multiply(2, q4.getVar(3), q4.getVar(4)));
		Polynomial<Fraction> p3 = q4.multiply(2, q4.getVar(1), q4.getVar(3));
		p3 = q4.subtract(p3, q4.multiply(2, q4.getVar(2), q4.getVar(4)));
		Polynomial<Fraction> p4 = q4.multiply(2, q4.getVar(1), q4.getVar(4));
		p4 = q4.add(p4, q4.multiply(2, q4.getVar(2), q4.getVar(3)));
		p4 = q4.subtract(p4, q4.multiply(2, q4.getVar(2), q4.getVar(4)));
		System.out.println(p1);
		System.out.println(p2);
		System.out.println(p3);
		System.out.println(p4);
		List<Polynomial<Fraction>> list = new ArrayList<>();
		list.add(p1);
		list.add(p2);
		list.add(p3);
		list.add(p4);
		List<AffinePoint<Fraction>> solution = q4.solve(list);
		System.out.println();
		for (AffinePoint<Fraction> point : solution) {
			System.out.println(point);
		}
		
		Ideal<Polynomial<Fraction>> idealOfInterest = q4.getIdeal(list);
		System.out.println();
		for (Polynomial<Fraction> polynomial : idealOfInterest.generators()) {
			System.out.println(polynomial);
		}
		System.exit(0);
		PAdicField fiveadic = new PAdicField(BigInteger.valueOf(5), 100);
		PAdicNumber five = fiveadic.getInteger(5);
		PAdicNumber four = fiveadic.getInteger(4);
		PAdicNumber minusone = fiveadic.getInteger(-1);
		PAdicNumber onefifth = fiveadic.inverse(fiveadic.getInteger(5));
		PAdicNumber onequarter = fiveadic.inverse(fiveadic.getInteger(4));
		PAdicNumber threequarter = fiveadic.divide(fiveadic.getInteger(3), fiveadic.getInteger(4));
		System.out.println("5 = " + five + " = " + fiveadic.toRational(five));
		System.out.println("4 = " + four + " = " + fiveadic.toRational(four));
		System.out.println("-1 = " + minusone + " = " + fiveadic.toRational(minusone));
		System.out.println("1/5 = " + onefifth + " = " + fiveadic.toRational(onefifth));
		System.out.println("1/4 = " + onequarter + " = " + fiveadic.toRational(onequarter));
		System.out.println("3/4 = " + threequarter + " = " + fiveadic.toRational(threequarter));

		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> qr = q.getUnivariatePolynomialRing();
		Polynomial<Fraction> poly = qr.getPolynomial(q.one(), q.zero(), q.zero(), q.zero(), q.zero(),
				q.zero(), q.zero(), q.one(), q.one());
		List<Polynomial<Fraction>> factors = q.factorization(poly);
		System.out.println("Factorization over the fractions:");
		for (Polynomial<Fraction> factor : factors) {
			System.out.println("factor: " + factor);
		}

		Reals r = Reals.r();
		Polynomial<Real> polyReal = MiscAlgorithms.mapPolynomial(poly, new MathMap<>() {

			@Override
			public Real evaluate(Fraction t) {
				return r.getEmbedding(t);
			}
		}, r.getUnivariatePolynomialRing());
		List<Polynomial<Real>> factorsReal = r.factorization(polyReal);
		System.out.println("Factorization over the reals:");
		for (Polynomial<Real> factor : factorsReal) {
			System.out.println("factor: " + factor);
		}

		Polynomial<Fraction> minpoly = qr.getPolynomial(q.one(), q.one(), q.one());
		NumberField nf = new NumberField(minpoly);
		Polynomial<NFE> polyNF = MiscAlgorithms.mapPolynomial(poly, new MathMap<>() {

			@Override
			public NFE evaluate(Fraction t) {
				return nf.getEmbedding(t);
			}
		}, nf.getUnivariatePolynomialRing());
		List<Polynomial<NFE>> factorsNF = nf.factorization(polyNF);
		System.out.println("Factorization over NumberField:");
		for (Polynomial<NFE> factor : factorsNF) {
			System.out.println("factor: " + factor);
		}

//		Map<Monomial, Fraction> intmap = new TreeMap<>(qr.getComparator());
//		intmap.put(qr.getMonomial(new int[] { 128 }), q.one());
//		intmap.put(qr.getMonomial(new int[] { 112 }), q.negative(q.one()));
//		intmap.put(qr.getMonomial(new int[] { 80 }), q.one());
//		intmap.put(qr.getMonomial(new int[] { 64 }), q.negative(q.one()));
//		intmap.put(qr.getMonomial(new int[] { 48 }), q.one());
//		intmap.put(qr.getMonomial(new int[] { 16 }), q.negative(q.one()));
//		intmap.put(qr.getMonomial(new int[] { 0 }), q.one());
//		Polynomial<Fraction> poly2 = qr.getPolynomial(intmap);
//		System.out.println("Factorization over Rational:");
//		List<Polynomial<Fraction>> factors2 = q.factorization(poly2);
//		for (Polynomial<Fraction> factor : factors2) {
//			System.out.println("factor: " + factor);
//		}

		ExtensionField<Fraction> extended = new ExtensionField<>(minpoly, q);
		MatrixAlgebra<Fraction> algebra = new FreeModule<>(q, 4).matrixAlgebra();
		List<List<Fraction>> m = new ArrayList<>();
		List<Fraction> row1 = new ArrayList<>();
		row1.add(q.getInteger(1));
		row1.add(q.getInteger(-1));
		row1.add(q.getInteger(11));
		row1.add(q.getInteger(-71));
		m.add(row1);
		List<Fraction> row2 = new ArrayList<>();
		row2.add(q.getInteger(1));
		row2.add(q.getInteger(-3));
		row2.add(q.getInteger(7));
		row2.add(q.getInteger(-49));
		m.add(row2);
		List<Fraction> row3 = new ArrayList<>();
		row3.add(q.getInteger(0));
		row3.add(q.getInteger(3));
		row3.add(q.getInteger(-10));
		row3.add(q.getInteger(35));
		m.add(row3);
		List<Fraction> row4 = new ArrayList<>();
		row4.add(q.getInteger(0));
		row4.add(q.getInteger(1));
		row4.add(q.getInteger(-10));
		row4.add(q.getInteger(37));
		m.add(row4);
		Matrix<Fraction> matrix = new Matrix<>(m);
		System.out.println(matrix);
		List<Fraction> rhs = new ArrayList<>();
		rhs.add(q.getInteger(1));
		rhs.add(q.getInteger(0));
		rhs.add(q.getInteger(0));
		rhs.add(q.getInteger(0));
		Vector<Fraction> vector = new Vector<>(rhs);
		System.out.println(vector);
		Vector<Fraction> result = algebra.solve(matrix, vector);
		System.out.println(result);
		ExtensionFieldElement<Fraction> alpha = extended.alpha();
		ExtensionFieldElement<Fraction> iE = extended.multiply(
				extended.getEmbedding(q.getFraction(z.getInteger(1), z.getInteger(24))),
				extended.add(extended.multiply(-127, alpha),
						extended.add(extended.multiply(-5, extended.power(alpha, 3)),
								extended.multiply(-19, extended.power(alpha, 5)),
								extended.multiply(-5, extended.power(alpha, 7)))));
		ExtensionFieldElement<Fraction> fourthRootOfTwo = extended.subtract(alpha, iE);
		ExtensionFieldElement<Fraction> sqrtOfTwo = extended.multiply(fourthRootOfTwo, fourthRootOfTwo);
		System.out.println(iE);
		System.out.println(fourthRootOfTwo);
		System.out.println(extended.multiply(iE, iE));
		System.out.println(extended.power(fourthRootOfTwo, 4));
		System.out.println("N((i-1)/sqrt(2)) = "
				+ extended.norm(extended.divide(extended.subtract(iE, extended.one()), sqrtOfTwo)));
		System.out.println("N(i) = " + extended.norm(iE));
		System.out.println("N(i-1) = " + extended.norm(extended.subtract(iE, extended.one())));
		System.out.println("N(4sqrt(2)) = " + extended.norm(fourthRootOfTwo));
		System.out.println("N(sqrt(2)) = " + extended.norm(sqrtOfTwo));
		System.out.println("N(2) = " + extended.norm(extended.getInteger(2)));
		System.out
				.println(
						"N((i-1-sqrt(2)-2) = "
								+ extended
										.norm(extended.divide(
												extended.subtract(iE,
														extended.add(extended.one(), sqrtOfTwo,
																extended.power(fourthRootOfTwo, 3))),
												extended.getInteger(2))));

		// System.exit(0);

		PrimeField base = new PrimeField(5);
		UnivariatePolynomialRing<PrimeFieldElement> ring = base.getUnivariatePolynomialRing();
		UnivariatePolynomial<PrimeFieldElement> mini = ring.getPolynomial(base.getElement(2), base.getElement(0),
				base.getElement(1));
		FiniteField f25 = new FiniteField(mini, base);
		FFE a = f25.getRandomElement();
		FFE b = f25.getRandomElement();
		// MatrixAlgebra<PrimeFieldElement> algebra = new MatrixAlgebra<>(base, 2);
		new TestMatrix<>(base, 2, 3);
		new TestMatrix<>(Integers.z(), 2, 3);
		new TestMatrix<>(Integers.z(), 3, 3);
		new TestMatrix<>(f25, 3, 3);
		new TestMatrix<>(f25.getUnivariatePolynomialRing(), 2, 3);
		new TestMatrix<>(base, 2, 2);
		new TestMatrix<>(Integers.z(), 2, 2);
		new TestMatrix<>(Integers.z(), 3, 2);
		new TestMatrix<>(f25, 3, 2);
		new TestMatrix<>(f25.getUnivariatePolynomialRing(), 2, 2);
		new TestEllipticCurve<>(f25, f25.one(), f25.one());
		new TestEllipticCurve<>(f25, f25.negative(f25.one()), f25.zero());

		mini = ring.getPolynomial(base.getElement(1), base.getElement(1), base.getElement(0), base.getElement(1));
		ExtensionField<PrimeFieldElement> f125 = new ExtensionField<PrimeFieldElement>(mini, base);
		Ideal<Polynomial<PrimeFieldElement>> ideal = ring.getIdeal(Collections.singletonList(mini));
		CoordinateRing<PrimeFieldElement> coord = new CoordinateRing<>(ring, ideal);
		EllipticCurve<ExtensionFieldElement<PrimeFieldElement>> ec = new EllipticCurve<>(f125, f125.one(), f125.one());

		System.out.println("Counting points!");
		System.out.println("Number of points: " + ec.getNumberOfElements());

		int n = 2;
		List<ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>>> torsion = ec.getTorsionPoints(n);
		System.out.println(torsion.size());
		for (ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point1 : torsion) {
			for (ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point2 : torsion) {
				System.out.println("P1: " + p1);
				System.out.println("P2: " + p2);
				BigInteger nBig = BigInteger.valueOf(n);
				System.out.println("Weil(P1, P2): " + ec.weilPairing(nBig, point1, point2));
				System.out.println("Weil(2P1, P2): " + ec.weilPairing(nBig, ec.multiply(2, point1), point2));
				System.out.println("Weil(P1, 2P2): " + ec.weilPairing(nBig, point1, ec.multiply(2, point2)));
			}
		}
		for (int i = 0; i < 20; i++) {
			ExtensionFieldElement<PrimeFieldElement> re = f125.getRandomElement();
			if (!re.equals(f125.zero())) {
				CoordinateRingElement<PrimeFieldElement> rp = coord.getEmbedding(f125.asPolynomial(re));
				System.out.println(re);
				System.out.println(rp);
				System.out.println(coord.isUnit(rp));
				System.out.println(f125.inverse(re));
				System.out.println(coord.inverse(rp));

				ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point1 = ec.getRandomElement();
				BigInteger order1 = ec.getOrder(point1);
				ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point2 = ec.getRandomElement();
				BigInteger order2 = ec.getOrder(point2);
				ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point3 = ec.getRandomElement();
				BigInteger order3 = ec.getOrder(point3);
				ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point4 = ec.getRandomElement();
				BigInteger order4 = ec.getOrder(point4);
				ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> point5 = ec.negative(ec.add(point2, point3));
				BigInteger order5 = ec.getOrder(point5);
				System.out.println("P1: " + point1 + " order: " + order1);
				System.out.println("P2: " + point2 + " order: " + order2);
				System.out.println("P3: " + point3 + " order: " + order3);
				System.out.println("P4: " + point4 + " order: " + order4);
				System.out.println("-(P2+P3): " + point5 + " order: " + order5);
				RationalFunction<ExtensionFieldElement<PrimeFieldElement>> f = ec.getRationalFunction(point1, point2,
						point3);
				System.out.println(f);
				f.simplify();
				System.out.println(f);
				FunctionField<ExtensionFieldElement<PrimeFieldElement>> ff = ec.getFunctionField();
				RationalFunction<ExtensionFieldElement<PrimeFieldElement>> f1 = ff.add(f, ff.one());
				RationalFunction<ExtensionFieldElement<PrimeFieldElement>> fsq = ff.multiply(f, f);
				System.out.println("f/f = " + ff.divide(f, f));
				System.out.println("f+f = " + ff.add(f, f));
				System.out.println("f*f = " + fsq);
				System.out.println("f+1 = " + f1);
				System.out.println("f(P1) = " + f.evaluate(point1));
				System.out.println("f(P2) = " + f.evaluate(point2));
				System.out.println("f(P3) = " + f.evaluate(point3));
				System.out.println("f(P4) = " + f.evaluate(point4));
				System.out.println("f(-(P2+P3)) = " + f.evaluate(point5));
				System.out.println("f(O) = " + f.evaluate(ec.neutral()));
				System.out.println("f+1(P1) = " + f1.evaluate(point1));
				System.out.println("f+1(P2) = " + f1.evaluate(point2));
				System.out.println("f+1(P3) = " + f1.evaluate(point3));
				System.out.println("f+1(P4) = " + f1.evaluate(point4));
				System.out.println("f+1(-(P2+P3)) = " + f1.evaluate(point5));
				System.out.println("f+1(O) = " + f1.evaluate(ec.neutral()));
				System.out.println("f*f(P1) = " + fsq.evaluate(point1));
				System.out.println("f*f(P2) = " + fsq.evaluate(point2));
				System.out.println("f*f(P3) = " + fsq.evaluate(point3));
				System.out.println("f*f(P4) = " + fsq.evaluate(point4));
				System.out.println("f*f(-(P2+P3)) = " + fsq.evaluate(point5));
				System.out.println("f*f(O) = " + fsq.evaluate(ec.neutral()));
			}
		}
		// new TestEllipticCurve<>(f125, f125.one(), f125.one());
		// new TestEllipticCurve<>(f125, f125.negative(f125.one()), f125.zero());
		/*
		 * TestAsymmetricCrypto test = new TestAsymmetricCrypto();
		 * test.doKeyAgreement(); test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 */
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger three = BigInteger.valueOf(3);
		PrimeField f54 = new PrimeField(two.pow(5).multiply(three.pow(4)).subtract(one));
		Polynomial<PrimeFieldElement> f54p = f54.getUnivariatePolynomialRing().getPolynomial(f54.one(), f54.zero(),
				f54.one());
		ExtensionField<PrimeFieldElement> f54sq = new ExtensionField<>(f54p, f54);
		new TestSike<ExtensionFieldElement<PrimeFieldElement>>(f54sq);

		PrimeField f9754 = new PrimeField(two.pow(97).multiply(three.pow(54)).subtract(one));
		Polynomial<PrimeFieldElement> f9754p = f9754.getUnivariatePolynomialRing().getPolynomial(f9754.one(),
				f9754.zero(), f9754.one());
		ExtensionField<PrimeFieldElement> f9754sq = new ExtensionField<>(f9754p, f9754);
		new TestSike<ExtensionFieldElement<PrimeFieldElement>>(f9754sq);

		PrimeField fp434 = new PrimeField(two.pow(0xd8).multiply(three.pow(0x89)).subtract(one));
		Polynomial<PrimeFieldElement> fp434p = fp434.getUnivariatePolynomialRing().getPolynomial(fp434.one(),
				fp434.zero(), fp434.one());
		ExtensionField<PrimeFieldElement> fp434sq = new ExtensionField<>(fp434p, fp434);
		new TestSike<ExtensionFieldElement<PrimeFieldElement>>(fp434sq);

		// BigPrimeField f65537 = new BigPrimeField(BigInteger.valueOf(65537));
		// ExtensionField<BigPrimeFieldElement> f65537sq = new ExtensionField<>(2,
		// f65537);
		// new TestSike<ExtensionFieldElement<BigPrimeFieldElement>>(f65537sq);
		// BigPrimeField f65537 = new BigPrimeField(BigInteger.valueOf(65537));
		// ExtensionField<BigPrimeFieldElement> f65537sq = new ExtensionField<>(2,
		// f65537);
		// new TestSike<ExtensionFieldElement<BigPrimeFieldElement>>(f65537sq);
		System.exit(0);
		// new TestEllipticCurve<>(f65537, f65537.one(), f65537.one());
		// new TestEllipticCurve<>(f65537, f65537.negative(f65537.one()),
		// f65537.zero());
		PrimeField f65537 = new PrimeField(BigInteger.valueOf(65537));
		TestFindRelatedCurves<PrimeFieldElement> f65537test = new TestFindRelatedCurves<>(f65537, f65537.one(),
				f65537.getInteger(2), new ProjectivePoint<PrimeFieldElement>(f65537, f65537.getInteger(1),
						f65537.getInteger(2), f65537.one()));
		f65537test.doAttack();

		// BigInteger p256Prime = new
		// BigInteger("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
		// 16);
		// BigInteger p256A = new BigInteger("-3");
		// BigInteger p256B = new BigInteger(
		// "41058363725152142129326129780047268409114441015993725554835256314039467401291");
		// BigPrimeField fp256 = new BigPrimeField(p256Prime);
		// BigInteger baseX = new BigInteger(
		// "48439561293906451759052585252797914202762949526041747995844080717082404635286");
		// BigInteger baseY = new BigInteger(
		// "36134250956749795798585127919587881956611106672985015071877198253568414405109");
		// ProjectivePoint<BigPrimeFieldElement> p256Generator = new
		// ProjectivePoint<>(fp256, fp256.getElement(baseX),
		// fp256.getElement(baseY), fp256.one());
		// new TestFindRelatedCurves(fp256, fp256.getElement(p256A),
		// fp256.getElement(p256B), p256Generator);
	}
}
