package main;

import java.math.BigInteger;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.CoordinateRing;
import fields.CoordinateRing.CoordinateRingElement;
import fields.ExtensionField;
import fields.ExtensionField.ExtensionFieldElement;
import fields.FunctionField;
import fields.IntegerRing;
import fields.MatrixAlgebra;
import fields.Polynomial;
import fields.Polynomial.Monomial;
import fields.PolynomialRing;
import fields.PrimeField;
import fields.PrimeField.PrimeFieldElement;
import fields.RationalFunction;
import varieties.ProjectivePoint;
import varieties.curves.EllipticCurve;

public class RunTests {

	public static void main(String[] args) {
		PrimeField base = new PrimeField(5);
		PolynomialRing<PrimeFieldElement> ring = ExtensionField.getPolynomialRing(base);
		Map<Monomial, PrimeFieldElement> map = new TreeMap<Monomial, PrimeFieldElement>(Polynomial.LEX);

		map.put(new Monomial(new int[] { 0 }), base.getElement(2));
		map.put(new Monomial(new int[] { 2 }), base.getElement(1));
		Polynomial<PrimeFieldElement> mini = new Polynomial<PrimeFieldElement>(map, ring);
		map.clear();
		ExtensionField<PrimeFieldElement> f25 = new ExtensionField<PrimeFieldElement>(mini, base);
		ExtensionFieldElement<PrimeFieldElement> a = f25.getRandomElement();
		ExtensionFieldElement<PrimeFieldElement> b = f25.getRandomElement();
		MatrixAlgebra<PrimeFieldElement> algebra = new MatrixAlgebra<>(base, 2);
		new TestMatrix<>(base, 2);
		new TestMatrix<>(new IntegerRing(), 2);
		new TestMatrix<>(new IntegerRing(), 3);
		new TestMatrix<>(f25, 3);
		new TestMatrix<>(new PolynomialRing<>(f25, 1, Polynomial.GREVLEX), 2);
		// new TestEllipticCurve<>(f25, f25.one(), f25.one());
		// new TestEllipticCurve<>(f25, f25.negative(f25.one()), f25.zero());

		map.put(new Monomial(new int[] { 0 }), base.getElement(1));
		map.put(new Monomial(new int[] { 1 }), base.getElement(1));
		map.put(new Monomial(new int[] { 3 }), base.getElement(1));
		mini = new Polynomial<PrimeFieldElement>(map, ring);
		ExtensionField<PrimeFieldElement> f125 = new ExtensionField<PrimeFieldElement>(mini, base);
		PolynomialRing<PrimeFieldElement>.Ideal ideal = ring.getIdeal(Collections.singletonList(mini));
		CoordinateRing<PrimeFieldElement> coord = new CoordinateRing<>(ring, ideal);
		EllipticCurve<ExtensionFieldElement<PrimeFieldElement>> ec = new EllipticCurve<>(f125, f125.one(), f125.one());

		System.out.println("Counting points!");
		System.out.println("Number of points: " + ec.getNumberOfElements());

		int n = 2;
		List<ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>>> torsion = ec.getTorsionPoints(n);
		System.out.println(torsion.size());
		for (ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> p1 : torsion) {
			for (ProjectivePoint<ExtensionFieldElement<PrimeFieldElement>> p2 : torsion) {
				System.out.println("P1: " + p1);
				System.out.println("P2: " + p2);
				BigInteger nBig = BigInteger.valueOf(n);
				System.out.println("Weil(P1, P2): " + ec.weilPairing(nBig, p1, p2));
				System.out.println("Weil(2P1, P2): " + ec.weilPairing(nBig, ec.multiply(2, p1), p2));
				System.out.println("Weil(P1, 2P2): " + ec.weilPairing(nBig, p1, ec.multiply(2, p2)));
			}
		}
		for (int i = 0; i < 20; i++) {
			ExtensionFieldElement<PrimeFieldElement> r = f125.getRandomElement();
			if (!r.equals(f125.zero())) {
				CoordinateRingElement<PrimeFieldElement> rp = coord.getEmbedding(f125.asPolynomial(r));
				System.out.println(r);
				System.out.println(rp);
				System.out.println(coord.isUnit(rp));
				System.out.println(f125.inverse(r));
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
		Map<Monomial, PrimeFieldElement> f54map = new TreeMap<>(Polynomial.LEX);
		f54map.put(new Monomial(new int[] { 0 }), f54.one());
		f54map.put(new Monomial(new int[] { 2 }), f54.one());
		Polynomial<PrimeFieldElement> f54p = new Polynomial<>(f54map, ExtensionField.getPolynomialRing(f54));
		ExtensionField<PrimeFieldElement> f54sq = new ExtensionField<>(f54p, f54);
		new TestSike<ExtensionFieldElement<PrimeFieldElement>>(f54sq);

		PrimeField f9754 = new PrimeField(two.pow(97).multiply(three.pow(54)).subtract(one));
		Map<Monomial, PrimeFieldElement> f9754map = new TreeMap<>(Polynomial.LEX);
		f9754map.put(new Monomial(new int[] { 0 }), f9754.one());
		f9754map.put(new Monomial(new int[] { 2 }), f9754.one());
		Polynomial<PrimeFieldElement> f9754p = new Polynomial<>(f9754map, ExtensionField.getPolynomialRing(f9754));
		ExtensionField<PrimeFieldElement> f9754sq = new ExtensionField<>(f9754p, f9754);
		new TestSike<ExtensionFieldElement<PrimeFieldElement>>(f9754sq);

		PrimeField fp434 = new PrimeField(two.pow(0xd8).multiply(three.pow(0x89)).subtract(one));
		Map<Monomial, PrimeFieldElement> fp434map = new TreeMap<>(Polynomial.LEX);
		fp434map.put(new Monomial(new int[] { 0 }), fp434.one());
		fp434map.put(new Monomial(new int[] { 2 }), fp434.one());
		Polynomial<PrimeFieldElement> fp434p = new Polynomial<>(fp434map, ExtensionField.getPolynomialRing(fp434));
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
